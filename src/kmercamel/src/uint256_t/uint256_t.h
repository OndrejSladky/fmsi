#pragma once
#ifndef __LITTLE_ENDIAN__
#ifndef __BIG_ENDIAN__
  #define __LITTLE_ENDIAN__ 1
#endif
#endif
#include "uint256_t.build"
#include <vector>
#include <cstring>

const uint128_t uint128_64(64);
const uint128_t uint128_128(128);
const uint128_t uint128_256(256);
const uint256_t uint256_0(0);
const uint256_t uint256_1(1);
const uint256_t uint256_max(uint128_t(-1), uint128_t(-1));

uint256_t::uint256_t(const bool & b)
    : uint256_t((uint8_t) b)
{}

uint256_t & uint256_t::operator=(const bool & rhs) {
    UPPER = 0;
    LOWER = rhs;
    return *this;
}

uint256_t::operator bool() const{
    return (bool) (UPPER | LOWER);
}

uint256_t::operator uint8_t() const{
    return (uint8_t) LOWER;
}

uint256_t::operator uint16_t() const{
    return (uint16_t) LOWER;
}

uint256_t::operator uint32_t() const{
    return (uint32_t) LOWER;
}

uint256_t::operator uint64_t() const{
    return (uint64_t) LOWER;
}

uint256_t::operator uint128_t() const{
    return LOWER;
}

uint256_t uint256_t::operator&(const uint128_t & rhs) const{
    return uint256_t(uint128_0, LOWER & rhs);
}

uint256_t uint256_t::operator&(const uint256_t & rhs) const{
    return uint256_t(UPPER & rhs.UPPER, LOWER & rhs.LOWER);
}

uint256_t & uint256_t::operator&=(const uint128_t & rhs){
    UPPER  = uint128_0;
    LOWER &= rhs;
    return *this;
}

uint256_t & uint256_t::operator&=(const uint256_t & rhs){
    UPPER &= rhs.UPPER;
    LOWER &= rhs.LOWER;
    return *this;
}

uint256_t uint256_t::operator|(const uint128_t & rhs) const{
    return uint256_t(UPPER , LOWER | rhs);
}

uint256_t uint256_t::operator|(const uint256_t & rhs) const{
    return uint256_t(UPPER | rhs.UPPER, LOWER | rhs.LOWER);
}

uint256_t & uint256_t::operator|=(const uint128_t & rhs){
    LOWER |= rhs;
    return *this;
}

uint256_t & uint256_t::operator|=(const uint256_t & rhs){
    UPPER |= rhs.UPPER;
    LOWER |= rhs.LOWER;
    return *this;
}

uint256_t uint256_t::operator^(const uint128_t & rhs) const{
    return uint256_t(UPPER, LOWER ^ rhs);
}

uint256_t uint256_t::operator^(const uint256_t & rhs) const{
    return uint256_t(UPPER ^ rhs.UPPER, LOWER ^ rhs.LOWER);
}

uint256_t & uint256_t::operator^=(const uint128_t & rhs){
    LOWER ^= rhs;
    return *this;
}

uint256_t & uint256_t::operator^=(const uint256_t & rhs){
    UPPER ^= rhs.UPPER;
    LOWER ^= rhs.LOWER;
    return *this;
}

uint256_t uint256_t::operator~() const{
    return uint256_t(~UPPER, ~LOWER);
}

uint256_t uint256_t::operator<<(const uint128_t & rhs) const{
    return *this << uint256_t(rhs);
}

uint256_t uint256_t::operator<<(const uint256_t & rhs) const{
    const uint128_t shift = rhs.LOWER;
    if (((bool) rhs.UPPER) || (shift >= uint128_256)){
        return uint256_0;
    }
    else if (shift == uint128_128){
        return uint256_t(LOWER, uint128_0);
    }
    else if (shift == uint128_0){
        return *this;
    }
    else if (shift < uint128_128){
        return uint256_t((UPPER << shift) + (LOWER >> (uint128_128 - shift)), LOWER << shift);
    }
    else if ((uint128_256 > shift) && (shift > uint128_128)){
        return uint256_t(LOWER << (shift - uint128_128), uint128_0);
    }
    else{
        return uint256_0;
    }
}

uint256_t & uint256_t::operator<<=(const uint128_t & shift){
    return *this <<= uint256_t(shift);
}

uint256_t & uint256_t::operator<<=(const uint256_t & shift){
    *this = *this << shift;
    return *this;
}

uint256_t uint256_t::operator>>(const uint128_t & rhs) const{
    return *this >> uint256_t(rhs);
}

uint256_t uint256_t::operator>>(const uint256_t & rhs) const{
    const uint128_t shift = rhs.LOWER;
    if (((bool) rhs.UPPER) | (shift >= uint128_256)){
        return uint256_0;
    }
    else if (shift == uint128_128){
        return uint256_t(UPPER);
    }
    else if (shift == uint128_0){
        return *this;
    }
    else if (shift < uint128_128){
        return uint256_t(UPPER >> shift, (UPPER << (uint128_128 - shift)) + (LOWER >> shift));
    }
    else if ((uint128_256 > shift) && (shift > uint128_128)){
        return uint256_t(UPPER >> (shift - uint128_128));
    }
    else{
        return uint256_0;
    }
}

uint256_t & uint256_t::operator>>=(const uint128_t & shift){
    return *this >>= uint256_t(shift);
}

uint256_t & uint256_t::operator>>=(const uint256_t & shift){
    *this = *this >> shift;
    return *this;
}

bool uint256_t::operator!() const{
    return ! (bool) *this;
}

bool uint256_t::operator&&(const uint128_t & rhs) const{
    return (*this && uint256_t(rhs));
}

bool uint256_t::operator&&(const uint256_t & rhs) const{
    return ((bool) *this && (bool) rhs);
}

bool uint256_t::operator||(const uint128_t & rhs) const{
    return (*this || uint256_t(rhs));
}

bool uint256_t::operator||(const uint256_t & rhs) const{
    return ((bool) *this || (bool) rhs);
}

bool uint256_t::operator==(const uint128_t & rhs) const{
    return (*this == uint256_t(rhs));
}

bool uint256_t::operator==(const uint256_t & rhs) const{
    return ((UPPER == rhs.UPPER) && (LOWER == rhs.LOWER));
}

bool uint256_t::operator!=(const uint128_t & rhs) const{
    return (*this != uint256_t(rhs));
}

bool uint256_t::operator!=(const uint256_t & rhs) const{
    return ((UPPER != rhs.UPPER) | (LOWER != rhs.LOWER));
}

bool uint256_t::operator>(const uint128_t & rhs) const{
    return (*this > uint256_t(rhs));
}

bool uint256_t::operator>(const uint256_t & rhs) const{
    if (UPPER == rhs.UPPER){
        return (LOWER > rhs.LOWER);
    }
    if (UPPER > rhs.UPPER){
        return true;
    }
    return false;
}

bool uint256_t::operator<(const uint128_t & rhs) const{
    return (*this < uint256_t(rhs));
}

bool uint256_t::operator<(const uint256_t & rhs) const{
    if (UPPER == rhs.UPPER){
        return (LOWER < rhs.LOWER);
    }
    if (UPPER < rhs.UPPER){
        return true;
    }
    return false;
}

bool uint256_t::operator>=(const uint128_t & rhs) const{
    return (*this >= uint256_t(rhs));
}

bool uint256_t::operator>=(const uint256_t & rhs) const{
    return ((*this > rhs) | (*this == rhs));
}

bool uint256_t::operator<=(const uint128_t & rhs) const{
    return (*this <= uint256_t(rhs));
}

bool uint256_t::operator<=(const uint256_t & rhs) const{
    return ((*this < rhs) | (*this == rhs));
}

uint256_t uint256_t::operator+(const uint128_t & rhs) const{
    return *this + uint256_t(rhs);
}

uint256_t uint256_t::operator+(const uint256_t & rhs) const{
    return uint256_t(UPPER + rhs.UPPER + (((LOWER + rhs.LOWER) < LOWER)?uint128_1:uint128_0), LOWER + rhs.LOWER);
}

uint256_t & uint256_t::operator+=(const uint128_t & rhs){
    return *this += uint256_t(rhs);
}

uint256_t & uint256_t::operator+=(const uint256_t & rhs){
    UPPER = rhs.UPPER + UPPER + ((LOWER + rhs.LOWER) < LOWER);
    LOWER = LOWER + rhs.LOWER;
    return *this;
}

uint256_t uint256_t::operator-(const uint128_t & rhs) const{
    return *this - uint256_t(rhs);
}

uint256_t uint256_t::operator-(const uint256_t & rhs) const{
    return uint256_t(UPPER - rhs.UPPER - ((LOWER - rhs.LOWER) > LOWER), LOWER - rhs.LOWER);
}

uint256_t & uint256_t::operator-=(const uint128_t & rhs){
    return *this -= uint256_t(rhs);
}

uint256_t & uint256_t::operator-=(const uint256_t & rhs){
    *this = *this - rhs;
    return *this;
}

uint256_t & uint256_t::operator++(){
    *this += uint256_1;
    return *this;
}

uint256_t uint256_t::operator++(int){
    uint256_t temp(*this);
    ++*this;
    return temp;
}

uint256_t & uint256_t::operator--(){
    *this -= uint256_1;
    return *this;
}

uint256_t uint256_t::operator--(int){
    uint256_t temp(*this);
    --*this;
    return temp;
}

uint256_t uint256_t::operator+() const{
    return *this;
}

uint256_t uint256_t::operator-() const{
    return ~*this + uint256_1;
}

const uint128_t & uint256_t::upper() const {
    return UPPER;
}

const uint128_t & uint256_t::lower() const {
    return LOWER;
}

uint256_t operator&(const uint128_t & lhs, const uint256_t & rhs){
    return rhs & lhs;
}

uint128_t & operator&=(uint128_t & lhs, const uint256_t & rhs){
    lhs = (rhs & lhs).lower();
    return lhs;
}

uint256_t operator|(const uint128_t & lhs, const uint256_t & rhs){
    return rhs | lhs;
}

uint128_t & operator|=(uint128_t & lhs, const uint256_t & rhs){
    lhs = (rhs | lhs).lower();
    return lhs;
}

uint256_t operator^(const uint128_t & lhs, const uint256_t & rhs){
    return rhs ^ lhs;
}

uint128_t & operator^=(uint128_t & lhs, const uint256_t & rhs){
    lhs = (rhs ^ lhs).lower();
    return lhs;
}

uint256_t operator<<(const bool & lhs, const uint256_t & rhs){
    return uint256_t(lhs) << rhs;
}

uint256_t operator<<(const uint8_t & lhs, const uint256_t & rhs){
    return uint256_t(lhs) << rhs;
}

uint256_t operator<<(const uint16_t & lhs, const uint256_t & rhs){
    return uint256_t(lhs) << rhs;
}

uint256_t operator<<(const uint32_t & lhs, const uint256_t & rhs){
    return uint256_t(lhs) << rhs;
}

uint256_t operator<<(const uint64_t & lhs, const uint256_t & rhs){
    return uint256_t(lhs) << rhs;
}

uint256_t operator<<(const uint128_t & lhs, const uint256_t & rhs){
    return uint256_t(lhs) << rhs;
}

uint256_t operator<<(const int8_t & lhs, const uint256_t & rhs){
    return uint256_t(lhs) << rhs;
}

uint256_t operator<<(const int16_t & lhs, const uint256_t & rhs){
    return uint256_t(lhs) << rhs;
}

uint256_t operator<<(const int32_t & lhs, const uint256_t & rhs){
    return uint256_t(lhs) << rhs;
}

uint256_t operator<<(const int64_t & lhs, const uint256_t & rhs){
    return uint256_t(lhs) << rhs;
}

uint128_t & operator<<=(uint128_t & lhs, const uint256_t & rhs){
    lhs = (uint256_t(lhs) << rhs).lower();
    return lhs;
}

uint256_t operator>>(const bool & lhs, const uint256_t & rhs){
    return uint256_t(lhs) >> rhs;
}

uint256_t operator>>(const uint8_t & lhs, const uint256_t & rhs){
    return uint256_t(lhs) >> rhs;
}

uint256_t operator>>(const uint16_t & lhs, const uint256_t & rhs){
    return uint256_t(lhs) >> rhs;
}

uint256_t operator>>(const uint32_t & lhs, const uint256_t & rhs){
    return uint256_t(lhs) >> rhs;
}

uint256_t operator>>(const uint64_t & lhs, const uint256_t & rhs){
    return uint256_t(lhs) >> rhs;
}

uint256_t operator>>(const uint128_t & lhs, const uint256_t & rhs){
    return uint256_t(lhs) >> rhs;
}

uint256_t operator>>(const int8_t & lhs, const uint256_t & rhs){
    return uint256_t(lhs) >> rhs;
}

uint256_t operator>>(const int16_t & lhs, const uint256_t & rhs){
    return uint256_t(lhs) >> rhs;
}

uint256_t operator>>(const int32_t & lhs, const uint256_t & rhs){
    return uint256_t(lhs) >> rhs;
}

uint256_t operator>>(const int64_t & lhs, const uint256_t & rhs){
    return uint256_t(lhs) >> rhs;
}

uint128_t & operator>>=(uint128_t & lhs, const uint256_t & rhs){
    lhs = (uint256_t(lhs) >> rhs).lower();
    return lhs;
}

// Comparison Operators
bool operator==(const uint128_t & lhs, const uint256_t & rhs){
    return rhs == lhs;
}

bool operator!=(const uint128_t & lhs, const uint256_t & rhs){
    return rhs != lhs;
}

bool operator>(const uint128_t & lhs, const uint256_t & rhs){
    return rhs < lhs;
}

bool operator<(const uint128_t & lhs, const uint256_t & rhs){
    return rhs > lhs;
}

bool operator>=(const uint128_t & lhs, const uint256_t & rhs){
    return rhs <= lhs;
}

bool operator<=(const uint128_t & lhs, const uint256_t & rhs){
    return rhs >= lhs;
}

// Arithmetic Operators
uint256_t operator+(const uint128_t & lhs, const uint256_t & rhs){
    return rhs + lhs;
}

uint128_t & operator+=(uint128_t & lhs, const uint256_t & rhs){
    lhs = (rhs + lhs).lower();
    return lhs;
}

uint256_t operator-(const uint128_t & lhs, const uint256_t & rhs){
    return -(rhs - lhs);
}

uint128_t & operator-=(uint128_t & lhs, const uint256_t & rhs){
    lhs = (-(rhs - lhs)).lower();
    return lhs;
}
