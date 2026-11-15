#pragma once

// The container of the files of FMSI 2 (format 5): one file per index, <name>.fmsi, and one per kLCP
// array, <name>.fmsi.klcp.
//
// A file is a prologue, sections, and the table of the sections at its end:
//
//   prologue (24 bytes)   magic "FMSI2IDX", u32 version (5), u32 number of sections, u64 offset of the table
//   sections              the bytes of each section, in any order
//   table                 one entry of 64 bytes per section (section_entry)
//
// All integers are little-endian. A section has a tag (up to 8 ASCII characters), an encoding (a number whose
// meaning belongs to the tag), flags, and its sizes. A reader looks up the sections it knows by their tags and
// skips the others, unless they are flagged as required: such a file needs a newer reader. New kinds of data
// are added as new tags, new representations of a kind as new encodings, so that the format grows without
// breaking the files that exist. The table is at the end, so that a section can be appended to a file.
//
// The section META is text, one "key=value" per line: what the file holds and how (see fms_index2.h). It is
// stored first and uncompressed, so that `head -c 1000 <file>` shows it.
//
// A section flagged DEFLATE is compressed in chunks of CHUNK_BYTES of its content, each a zlib stream of its
// own, so that the chunks are compressed and decompressed on several threads and a reader never holds more than
// a few chunks: the chunks, then their compressed sizes (u32 each), the number of chunks (u64) and the size
// of a chunk (u32).

#include <atomic>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <functional>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <zlib.h>

#include "parallel.h"

namespace fmsi2 {
namespace file {

static const char MAGIC[8] = {'F', 'M', 'S', 'I', '2', 'I', 'D', 'X'};
static const uint32_t CONTAINER_VERSION = 5;
/// Size of a chunk of a compressed section and of a piece read at a time: small, so that loading an index holds
/// little beside the index (a reader on t threads holds t chunks, compressed and not).
static const uint32_t CHUNK_BYTES = 1 << 18;

enum : uint32_t {
    /// A reader that does not know the tag or the encoding of the section must refuse the file.
    REQUIRED = 1,
    /// The section is compressed in chunks by zlib.
    DEFLATE = 2,
};

struct section_entry {
    char tag[8];
    uint32_t encoding;
    uint32_t flags;
    /// Position of the section in the file and its size there.
    uint64_t offset;
    uint64_t stored_bytes;
    /// Size of the content (after decompression).
    uint64_t raw_bytes;
    /// Number of items of the content (characters, bits), whose meaning belongs to the tag.
    uint64_t items;
    /// CRC-32 of the content.
    uint32_t crc32;
    uint32_t reserved32;
    uint64_t reserved64;

    std::string name() const { return std::string(tag, strnlen(tag, sizeof(tag))); }
};
static_assert(sizeof(section_entry) == 64, "a table entry has 64 bytes");

inline void require_little_endian() {
    const uint32_t one = 1;
    if (*reinterpret_cast<const uint8_t *>(&one) != 1) {
        throw std::runtime_error("the files of FMSI 2 are little-endian; this machine is not");
    }
}

/// Key-value metadata of a file, kept in the order of insertion.
class metadata {
    std::vector<std::pair<std::string, std::string>> items_;

  public:
    void set(const std::string &key, const std::string &value) {
        for (auto &item : items_) {
            if (item.first == key) {
                item.second = value;
                return;
            }
        }
        items_.emplace_back(key, value);
    }

    void set(const std::string &key, uint64_t value) { set(key, std::to_string(value)); }

    bool has(const std::string &key) const {
        for (const auto &item : items_) {
            if (item.first == key) return true;
        }
        return false;
    }

    std::string get(const std::string &key, const std::string &fallback = "") const {
        for (const auto &item : items_) {
            if (item.first == key) return item.second;
        }
        return fallback;
    }

    /// The value as a number; it must be there.
    uint64_t number(const std::string &key) const {
        if (!has(key)) throw std::invalid_argument("the file has no '" + key + "' in its metadata");
        return std::stoull(get(key));
    }

    uint64_t number(const std::string &key, uint64_t fallback) const { return has(key) ? std::stoull(get(key)) : fallback; }

    const std::vector<std::pair<std::string, std::string>> &items() const { return items_; }

    std::string text() const {
        std::string result;
        for (const auto &item : items_) result += item.first + "=" + item.second + "\n";
        return result;
    }

    static metadata parse(const std::string &text) {
        metadata result;
        std::istringstream lines(text);
        std::string line;
        while (std::getline(lines, line)) {
            size_t eq = line.find('=');
            if (eq == std::string::npos) continue;
            result.set(line.substr(0, eq), line.substr(eq + 1));
        }
        return result;
    }
};

class writer {
    std::ofstream out_;
    std::string path_;
    std::vector<section_entry> table_;

    static section_entry entry(const std::string &tag, uint32_t encoding, uint32_t flags, uint64_t items) {
        if (tag.empty() || tag.size() > 8) throw std::invalid_argument("a section tag has 1 to 8 characters");
        section_entry e;
        memset(&e, 0, sizeof(e));
        memcpy(e.tag, tag.data(), tag.size());
        e.encoding = encoding;
        e.flags = flags;
        e.items = items;
        return e;
    }

  public:
    explicit writer(const std::string &path) : out_(path, std::ios::binary), path_(path) {
        require_little_endian();
        if (!out_) throw std::invalid_argument("couldn't open file " + path);
        // The prologue is completed by finish().
        char prologue[24] = {0};
        out_.write(prologue, sizeof(prologue));
    }

    /// A section stored as it is.
    void add(const std::string &tag, uint32_t encoding, uint32_t flags, uint64_t items, const void *data, size_t bytes) {
        section_entry e = entry(tag, encoding, flags & ~(uint32_t)DEFLATE, items);
        e.offset = (uint64_t)out_.tellp();
        e.stored_bytes = e.raw_bytes = bytes;
        e.crc32 = (uint32_t)crc32_z(0, reinterpret_cast<const Bytef *>(data), bytes);
        out_.write(reinterpret_cast<const char *>(data), bytes);
        table_.push_back(e);
    }

    /// A section of [bytes] bytes produced piece by piece: fill(begin, end, buffer) writes the bytes [begin, end)
    /// of the content to the buffer; it is called from several threads, for consecutive pieces of CHUNK_BYTES.
    /// With [deflate], the pieces are compressed.
    void add_stream(const std::string &tag, uint32_t encoding, uint32_t flags, uint64_t items, uint64_t bytes, bool deflate,
                    unsigned threads, const std::function<void(uint64_t, uint64_t, uint8_t *)> &fill, int level = 6) {
        section_entry e = entry(tag, encoding, deflate ? (flags | DEFLATE) : (flags & ~(uint32_t)DEFLATE), items);
        e.offset = (uint64_t)out_.tellp();
        e.raw_bytes = bytes;
        const uint64_t chunks = (bytes + CHUNK_BYTES - 1) / CHUNK_BYTES;
        std::vector<uint32_t> sizes(chunks);
        std::vector<uint32_t> crcs(chunks);
        // A batch of chunks is produced and compressed on the threads, then written in its order.
        const uint64_t batch = std::max<uint64_t>(1, threads) * 4;
        std::vector<std::vector<uint8_t>> raw(batch), packed(batch);
        for (uint64_t first = 0; first < chunks; first += batch) {
            const uint64_t count = std::min(batch, chunks - first);
            std::atomic<bool> failed{false};
            parallel_chunks(count, threads, 1, [&](size_t begin, size_t end, unsigned) {
                for (size_t c = begin; c < end; ++c) {
                    const uint64_t from = (first + c) * CHUNK_BYTES, to = std::min(bytes, from + CHUNK_BYTES);
                    raw[c].resize(to - from);
                    fill(from, to, raw[c].data());
                    crcs[first + c] = (uint32_t)crc32_z(0, raw[c].data(), raw[c].size());
                    if (deflate) {
                        uLongf size = compressBound(raw[c].size());
                        packed[c].resize(size);
                        if (compress2(packed[c].data(), &size, raw[c].data(), raw[c].size(), level) != Z_OK) failed = true;
                        packed[c].resize(size);
                    }
                }
            });
            if (failed) throw std::runtime_error("zlib failed to compress a chunk of " + path_);
            for (uint64_t c = 0; c < count; ++c) {
                const std::vector<uint8_t> &data = deflate ? packed[c] : raw[c];
                sizes[first + c] = (uint32_t)data.size();
                out_.write(reinterpret_cast<const char *>(data.data()), data.size());
                e.crc32 = (uint32_t)crc32_combine(e.crc32, crcs[first + c], raw[c].size());
            }
        }
        if (deflate) {
            out_.write(reinterpret_cast<const char *>(sizes.data()), sizes.size() * sizeof(uint32_t));
            out_.write(reinterpret_cast<const char *>(&chunks), sizeof(chunks));
            out_.write(reinterpret_cast<const char *>(&CHUNK_BYTES), sizeof(CHUNK_BYTES));
        }
        e.stored_bytes = (uint64_t)out_.tellp() - e.offset;
        table_.push_back(e);
    }

    /// Write the table and the prologue; the file is complete only after this.
    void finish() {
        const uint64_t table_offset = (uint64_t)out_.tellp();
        out_.write(reinterpret_cast<const char *>(table_.data()), table_.size() * sizeof(section_entry));
        out_.seekp(0);
        out_.write(MAGIC, sizeof(MAGIC));
        const uint32_t version = CONTAINER_VERSION, count = (uint32_t)table_.size();
        out_.write(reinterpret_cast<const char *>(&version), sizeof(version));
        out_.write(reinterpret_cast<const char *>(&count), sizeof(count));
        out_.write(reinterpret_cast<const char *>(&table_offset), sizeof(table_offset));
        out_.close();
        if (!out_) throw std::invalid_argument("couldn't write file " + path_);
    }

    const std::vector<section_entry> &sections() const { return table_; }
};

/// Whether the file starts with the magic of FMSI 2, and its version (0 if it does not).
inline uint32_t file_version(const std::string &path) {
    std::ifstream in(path, std::ios::binary);
    char magic[sizeof(MAGIC)];
    uint32_t version = 0;
    in.read(magic, sizeof(magic));
    in.read(reinterpret_cast<char *>(&version), sizeof(version));
    if (!in || memcmp(magic, MAGIC, sizeof(MAGIC)) != 0) return 0;
    return version;
}

class reader {
    mutable std::ifstream in_;
    std::string path_;
    std::vector<section_entry> table_;
    uint64_t file_bytes_ = 0;

  public:
    explicit reader(const std::string &path) : in_(path, std::ios::binary), path_(path) {
        require_little_endian();
        if (!in_) throw std::invalid_argument("couldn't open file " + path);
        char magic[sizeof(MAGIC)];
        uint32_t version = 0, count = 0;
        uint64_t table_offset = 0;
        in_.read(magic, sizeof(magic));
        in_.read(reinterpret_cast<char *>(&version), sizeof(version));
        in_.read(reinterpret_cast<char *>(&count), sizeof(count));
        in_.read(reinterpret_cast<char *>(&table_offset), sizeof(table_offset));
        if (!in_ || memcmp(magic, MAGIC, sizeof(MAGIC)) != 0) {
            throw std::invalid_argument("file " + path + " is not a file of FMSI 2");
        }
        if (version != CONTAINER_VERSION) {
            throw std::invalid_argument("file " + path + " has format version " + std::to_string(version) +
                                        ", this reader knows version " + std::to_string(CONTAINER_VERSION));
        }
        in_.seekg(0, std::ios::end);
        file_bytes_ = (uint64_t)in_.tellg();
        if (table_offset > file_bytes_ || (file_bytes_ - table_offset) / sizeof(section_entry) < count) {
            throw std::invalid_argument("file " + path + " is truncated");
        }
        table_.resize(count);
        in_.seekg((std::streamoff)table_offset);
        in_.read(reinterpret_cast<char *>(table_.data()), count * sizeof(section_entry));
        if (!in_) throw std::invalid_argument("file " + path + " is truncated");
        for (const auto &e : table_) {
            if (e.offset > file_bytes_ || e.stored_bytes > file_bytes_ - e.offset) {
                throw std::invalid_argument("file " + path + " is truncated");
            }
        }
    }

    const std::string &path() const { return path_; }

    uint64_t file_bytes() const { return file_bytes_; }

    const std::vector<section_entry> &sections() const { return table_; }

    const section_entry *find(const std::string &tag) const {
        for (const auto &e : table_) {
            if (e.name() == tag) return &e;
        }
        return nullptr;
    }

    const section_entry &get(const std::string &tag) const {
        const section_entry *e = find(tag);
        if (!e) throw std::invalid_argument("file " + path_ + " has no section " + tag);
        return *e;
    }

    /// Refuse the file if it has a required section whose tag is not among the known ones.
    void check_known(const std::vector<std::string> &known) const {
        for (const auto &e : table_) {
            if (!(e.flags & REQUIRED)) continue;
            bool found = false;
            for (const auto &tag : known) found = found || tag == e.name();
            if (!found) {
                throw std::invalid_argument("file " + path_ + " has a section (" + e.name() +
                                            ") that this version of FMSI does not know; a newer one is needed");
            }
        }
    }

    /// The content of a section, pieces of it in their order: consume(offset, data, bytes), called from one
    /// thread. The chunks of a compressed section are decompressed on the given threads.
    void read(const section_entry &e, unsigned threads, const std::function<void(uint64_t, const uint8_t *, size_t)> &consume,
              bool verify = false) const {
        uint32_t crc = 0;
        if (!(e.flags & DEFLATE)) {
            std::vector<uint8_t> buffer(std::min<uint64_t>(e.stored_bytes, (uint64_t)CHUNK_BYTES));
            in_.seekg((std::streamoff)e.offset);
            for (uint64_t done = 0; done < e.stored_bytes;) {
                const size_t bytes = (size_t)std::min<uint64_t>(buffer.size(), e.stored_bytes - done);
                in_.read(reinterpret_cast<char *>(buffer.data()), bytes);
                if (!in_) throw std::invalid_argument("file " + path_ + " is truncated");
                if (verify) crc = (uint32_t)crc32_z(crc, buffer.data(), bytes);
                consume(done, buffer.data(), bytes);
                done += bytes;
            }
        } else {
            uint64_t chunks = 0;
            uint32_t chunk_bytes = 0;
            if (e.stored_bytes < sizeof(chunks) + sizeof(chunk_bytes)) throw std::invalid_argument("file " + path_ + " is damaged");
            in_.seekg((std::streamoff)(e.offset + e.stored_bytes - sizeof(chunks) - sizeof(chunk_bytes)));
            in_.read(reinterpret_cast<char *>(&chunks), sizeof(chunks));
            in_.read(reinterpret_cast<char *>(&chunk_bytes), sizeof(chunk_bytes));
            if (!in_ || chunk_bytes == 0 || chunks != (e.raw_bytes + chunk_bytes - 1) / chunk_bytes ||
                chunks > (e.stored_bytes - sizeof(chunks) - sizeof(chunk_bytes)) / sizeof(uint32_t)) {
                throw std::invalid_argument("file " + path_ + " is damaged");
            }
            std::vector<uint32_t> sizes(chunks);
            in_.seekg((std::streamoff)(e.offset + e.stored_bytes - sizeof(chunks) - sizeof(chunk_bytes) - chunks * sizeof(uint32_t)));
            in_.read(reinterpret_cast<char *>(sizes.data()), chunks * sizeof(uint32_t));
            in_.seekg((std::streamoff)e.offset);
            const uint64_t batch = std::max<uint64_t>(1, threads);
            std::vector<std::vector<uint8_t>> packed(batch), raw(batch);
            for (uint64_t first = 0; first < chunks; first += batch) {
                const uint64_t count = std::min(batch, chunks - first);
                for (uint64_t c = 0; c < count; ++c) {
                    packed[c].resize(sizes[first + c]);
                    in_.read(reinterpret_cast<char *>(packed[c].data()), packed[c].size());
                }
                if (!in_) throw std::invalid_argument("file " + path_ + " is truncated");
                // A thread must not throw: the failure is reported after the batch.
                std::atomic<bool> damaged{false};
                parallel_chunks(count, threads, 1, [&](size_t begin, size_t end, unsigned) {
                    for (size_t c = begin; c < end; ++c) {
                        const uint64_t from = (first + c) * chunk_bytes;
                        raw[c].resize((size_t)std::min<uint64_t>(chunk_bytes, e.raw_bytes - from));
                        uLongf size = raw[c].size();
                        if (uncompress(raw[c].data(), &size, packed[c].data(), packed[c].size()) != Z_OK || size != raw[c].size()) {
                            damaged = true;
                        }
                    }
                });
                if (damaged) throw std::invalid_argument("file " + path_ + " is damaged (section " + e.name() + ")");
                for (uint64_t c = 0; c < count; ++c) {
                    if (verify) crc = (uint32_t)crc32_z(crc, raw[c].data(), raw[c].size());
                    consume((first + c) * chunk_bytes, raw[c].data(), raw[c].size());
                }
            }
        }
        if (verify && crc != e.crc32) {
            throw std::invalid_argument("file " + path_ + " is damaged: the checksum of section " + e.name() + " does not match");
        }
    }

    /// The whole content of a (small) section.
    std::string read_all(const section_entry &e, bool verify = true) const {
        std::string content((size_t)e.raw_bytes, '\0');
        read(e, 1, [&](uint64_t offset, const uint8_t *data, size_t bytes) { memcpy(&content[(size_t)offset], data, bytes); }, verify);
        return content;
    }

    metadata meta() const { return metadata::parse(read_all(get("META"))); }
};

} // namespace file
} // namespace fmsi2
