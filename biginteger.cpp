#include<iostream>
#include"biginteger.h"
#include<string>
constexpr long long two_in_32 = 1ll<<32;
constexpr int additional_base = 1000000000;

void delete_zeroes_str(std::string& number) {
	while (number[0] == '0' && number.size() > 1) {
		number.erase(number.begin(), number.begin() + 1);
	}
}
void BigInteger::delete_zeroes() {
	unsigned int counter_zeroes = 0;
	while (m_number_of_bits_int32 - counter_zeroes > 1 && m_bits_of_int32[m_number_of_bits_int32 - 1 - counter_zeroes] == 0) {
		++counter_zeroes;
	}
	if (counter_zeroes) {
		m_number_of_bits_int32 -= counter_zeroes;
		unsigned int* copy = new unsigned int[m_number_of_bits_int32];
		for (size_t i = 0; i < m_number_of_bits_int32; ++i) {
			copy[i] = m_bits_of_int32[i];
		}
		delete[] m_bits_of_int32;
		m_bits_of_int32 = copy;
		copy = nullptr;
		delete[] copy;
	}
}

BigInteger::BigInteger() : m_sign(true), m_number_of_bits_int32(1) {
	m_bits_of_int32 = new unsigned int[1];
	m_bits_of_int32[0] = 0;
}
BigInteger::BigInteger(long long number) {
	if (number < 0) {
		m_sign = false;
		number *= -1;
	}
	else {
		m_sign = true;
	}
	if (number < two_in_32) {
		m_bits_of_int32 = new unsigned int[1];
		m_bits_of_int32[0] =  number;
		m_number_of_bits_int32 = 1;
	}
	else {
		m_number_of_bits_int32 = 2;
		m_bits_of_int32 = new unsigned int[2];
		m_bits_of_int32[1] = number / two_in_32;
		m_bits_of_int32[0] = number % two_in_32;

	}
}
BigInteger::BigInteger(unsigned long long number) {
		m_sign = true;
	if (number < two_in_32) {
		m_bits_of_int32 = new unsigned int[1];
		m_bits_of_int32[0] = number;
		m_number_of_bits_int32 = 1;
	}
	else {
		m_number_of_bits_int32 = 2;
		m_bits_of_int32 = new unsigned int[2];
		m_bits_of_int32[1] = number / two_in_32;
		m_bits_of_int32[0] = number % two_in_32;

	}
}
BigInteger::BigInteger(const BigInteger& X) : m_sign(X.m_sign), m_number_of_bits_int32(X.m_number_of_bits_int32) {
	m_bits_of_int32 = new unsigned int[X.m_number_of_bits_int32];
	for (size_t i = 0; i < m_number_of_bits_int32; ++i) 
		m_bits_of_int32[i] = X.m_bits_of_int32[i];
}
BigInteger::BigInteger(BigInteger&& X) noexcept {
	m_sign = X.m_sign;
	m_number_of_bits_int32 = X.m_number_of_bits_int32;
	m_bits_of_int32 = X.m_bits_of_int32;
	X.m_sign = true;
	X.m_number_of_bits_int32 = 0;
	X.m_bits_of_int32 = nullptr;
}
BigInteger::BigInteger(std::string number) {
	if (number.size() == 0) {
		throw std::invalid_argument(number);
	}
	else if (number[0] == '-') {
		this->m_sign = false;
		number.erase(number.begin(), number.begin() + 1);
		delete_zeroes_str(number);
		if (number.size() == 0)
			throw std::invalid_argument(number);
	}
	else {
		delete_zeroes_str(number);
		this->m_sign = true;
	}
	if (!CheckString(number))
		throw std::invalid_argument(number);
	m_number_of_bits_int32 = 0;
	m_bits_of_int32 = new unsigned int[0];
	processing_string(number);
}

void BigInteger::processing_string(std::string& number) {
	unsigned int size_new=0;
	if (number.length() % 9 == 0) {
		size_new = number.length()/9;
	}
	else {
		size_new = number.length() / 9 + 1;
	}
	unsigned long long int* bits_new_base = new unsigned long long[size_new];
	for (int i = number.length(), j=0; i > 0; i -= 9, ++j) {
		if (i < 9) {
			bits_new_base[j] = atoi(number.substr(0, i).c_str());
		}
		else {
			bits_new_base[j] = atoi(number.substr(i - 9, 9).c_str());
		}
	}
	unsigned int counter_zeroes = 0;
	while (m_number_of_bits_int32 - counter_zeroes > 1 && m_bits_of_int32[m_number_of_bits_int32 - 1 - counter_zeroes] == 0) {
		++counter_zeroes;
	}
	if (counter_zeroes) {
		size_new -= counter_zeroes;
		unsigned long long* copy = new unsigned long long[size_new];
		for (size_t i = 0; i < size_new; ++i) {
			copy[i] = bits_new_base[i];
		}
		delete[] bits_new_base;
		bits_new_base = copy;
		copy = nullptr;
		delete[] copy;
	}
	do {
		unsigned long long carry = 0;
		for (int i = 0; i < size_new; ++i) {
			unsigned long long cur = carry * additional_base + bits_new_base[size_new - i - 1];
			bits_new_base[size_new - i - 1] = cur / two_in_32;
			carry = cur % two_in_32;
		}
		this->resize_up();
		this->m_bits_of_int32[m_number_of_bits_int32 - 1] = carry;
	} while (!CheckZero(size_new, bits_new_base));
}

bool CheckZero(int size_new, unsigned long long* bits_new_base) {
	for (int i = 0; i <size_new; ++i) {
		if (*(bits_new_base+i) != 0)
			return false;
	}
	return true;
}

BigInteger::~BigInteger() {
	delete[] m_bits_of_int32;
}

BigInteger& BigInteger::operator=(const BigInteger& X) {
	if (this != &X) {
		m_sign = X.m_sign;
		m_number_of_bits_int32 = X.m_number_of_bits_int32;
		delete[] m_bits_of_int32;
		m_bits_of_int32 = new unsigned int[X.m_number_of_bits_int32];
		for (size_t i = 0; i < m_number_of_bits_int32; ++i)
			m_bits_of_int32[i] = X.m_bits_of_int32[i];
	}
	return *this;
}

BigInteger& BigInteger::operator=(BigInteger&& X) noexcept {
	if (this != &X) {
		delete[] m_bits_of_int32;
		m_sign = X.m_sign;
		m_number_of_bits_int32 = X.m_number_of_bits_int32;
		m_bits_of_int32 = X.m_bits_of_int32;
		X.m_bits_of_int32 = nullptr;
		X.m_number_of_bits_int32 = 0;
	}
	return *this;
}

bool BigInteger::get_sign() const {
	return m_sign;
}

unsigned int BigInteger::get_number_of_bits_int32() const {
	return m_number_of_bits_int32;
}

unsigned int BigInteger::operator[](int pos) const {
	return this->m_bits_of_int32[pos];
}

bool operator==(const BigInteger& first, const BigInteger& second) {
	if (!first.m_sign && second.m_sign || first.m_sign && !second.m_sign || first.m_number_of_bits_int32!= second.m_number_of_bits_int32)
		return false;
	for (size_t i = 0; i < first.m_number_of_bits_int32; ++i) {
		if (first.m_bits_of_int32[i] != second.m_bits_of_int32[i]) 
			return false;
	}
	return true;
}

bool operator!=(const BigInteger& first, const BigInteger& second) {
	return !(first==second);
}

bool operator<(const BigInteger& first, const BigInteger& second){
	if (first.m_sign && !second.m_sign) return false;
	if (!first.m_sign && second.m_sign) return true;
	if (first.m_sign) {
		if (first.m_number_of_bits_int32 < second.m_number_of_bits_int32) return true;
		if (first.m_number_of_bits_int32 > second.m_number_of_bits_int32) return false;
		for (int i = first.m_number_of_bits_int32-1; i > -1; --i) {
			if (first.m_bits_of_int32[i] != second.m_bits_of_int32[i])
				return first.m_bits_of_int32[i] < second.m_bits_of_int32[i];
		}
		return false;
	}
	else {
		if (first.m_number_of_bits_int32 > second.m_number_of_bits_int32) return true;
		if (first.m_number_of_bits_int32 < second.m_number_of_bits_int32) return false;
		for (int i = first.m_number_of_bits_int32-1; i > -1; --i) {
			if (first.m_bits_of_int32[i] != second.m_bits_of_int32[i]) 
				return first.m_bits_of_int32[i] > second.m_bits_of_int32[i];
		}
		return false;
	}
}

bool operator>(const BigInteger& first, const BigInteger& second) {
	return second < first;
}

bool operator>=(const BigInteger& first, const BigInteger& second) {
	return !(first < second);
}

bool operator<=(const BigInteger& first, const BigInteger& second) {
	return (first < second || first==second);
}

bool cmp_module(const BigInteger& first, const BigInteger& second) {
	if (first.m_number_of_bits_int32 > second.m_number_of_bits_int32)
		return true;
	if (first.m_number_of_bits_int32 < second.m_number_of_bits_int32)
		return false;
	unsigned int counter_eq_bits = 0;
	for (int i = first.m_number_of_bits_int32; i > -1; --i) {
		if (first.m_bits_of_int32[i] > second.m_bits_of_int32[i]) return true;
		if (first.m_bits_of_int32[i] == second.m_bits_of_int32[i]) ++counter_eq_bits;
	}
	if (++counter_eq_bits == first.m_number_of_bits_int32) return true;
	return false;
}

const BigInteger BigInteger::operator+() const {
	return *this;
}

const BigInteger BigInteger::operator-() const{
	BigInteger copy = *this;
	if (copy.m_bits_of_int32[0] == 0 && copy.m_number_of_bits_int32 == 1)
		return copy;
	copy.m_sign = !copy.m_sign;
	return copy;
}

BigInteger& BigInteger::operator+=(const BigInteger& second) {
	if (this->m_sign && second.m_sign || !this->m_sign && !second.m_sign) {
		*this = std::move(add_numbers(*this,second));
	}
	else if (cmp_module(*this, second) && (this->m_sign && !second.m_sign || !this->m_sign && second.m_sign)) {
		substract_numbers(*this,second);
	}
	else if (!cmp_module(*this, second) && (this->m_sign && !second.m_sign || !this->m_sign && second.m_sign)) {
		BigInteger copy_1 = second, copy_2=*this;
		*this=substract_numbers(copy_1, copy_2);
		this->m_sign != this->m_sign;
	}
	return *this;
}

BigInteger& BigInteger::operator-=(const BigInteger& second) {
	if (*this >= second && this->m_sign && second.m_sign) {
		substract_numbers(*this, second);
	}
	else if (second >= *this && this->m_sign && second.m_sign || *this >= second && !this->m_sign && !second.m_sign) {
		BigInteger copy_1 = second, copy_2=*this;
		*this = substract_numbers(copy_1, copy_2);
		this->m_sign = !this->m_sign;
	}
	else if (this->m_sign && !second.m_sign || !this->m_sign && second.m_sign) {
		add_numbers(*this,second);
	}
	else if (*this < second && !this->m_sign && !second.m_sign) {
		substract_numbers(*this, second);
	}
	return *this;
}

void BigInteger::resize(unsigned int size) {
    delete[] m_bits_of_int32;
	m_bits_of_int32 = new unsigned int[size];
    m_number_of_bits_int32 = size;
	for (size_t i = 0; i < size; ++i)
		m_bits_of_int32[i] = 0;
}

void BigInteger::resize_down() {
	--this->m_number_of_bits_int32;
	unsigned int* copy_bits_of_int32 = new unsigned int[this->m_number_of_bits_int32];
	for (size_t i = 0; i < this->m_number_of_bits_int32; ++i) {
		copy_bits_of_int32[i] = this->m_bits_of_int32[i];
	}
	delete[] m_bits_of_int32;
	m_bits_of_int32 = copy_bits_of_int32;
	copy_bits_of_int32 = nullptr;
	delete[] copy_bits_of_int32;
}
void BigInteger::resize_up() {
	++this->m_number_of_bits_int32;
	unsigned int* copy_bits_of_int32 = new unsigned int[this->m_number_of_bits_int32];
	for (size_t i = 0; i < this->m_number_of_bits_int32-1; ++i) {
		copy_bits_of_int32[i] = this->m_bits_of_int32[i];
	}
	copy_bits_of_int32[this->m_number_of_bits_int32-1] = 0;
	delete[] m_bits_of_int32;
	m_bits_of_int32 = copy_bits_of_int32;
	copy_bits_of_int32 = nullptr;
	delete[] copy_bits_of_int32;
}

BigInteger& add_numbers(BigInteger& first, const BigInteger& second) {
	BigInteger tmp(0);
	tmp.m_sign = first.m_sign;
	unsigned int size_tmp = (first.m_number_of_bits_int32 > second.m_number_of_bits_int32 ? 
		first.m_number_of_bits_int32 + 1 : second.m_number_of_bits_int32+1);
	tmp.resize(size_tmp);
	unsigned long long sum = 0;
	int carry = 0;
	size_t i;
	for (i = 0; i < size_tmp-1 || carry; ++i) {
		sum = (unsigned long long)carry + (i < first.m_number_of_bits_int32 ? first.m_bits_of_int32[i] : 0)
			+ (i < second.m_number_of_bits_int32 ? second.m_bits_of_int32[i] : 0);
		carry = (sum  >= two_in_32);
		if (carry) tmp.m_bits_of_int32[i] = sum - two_in_32;
		else tmp.m_bits_of_int32[i] = sum;
	}
    if(i != size_tmp){
        tmp.resize_down();
    }
	first = tmp;
	return first;
}

BigInteger& substract_numbers(BigInteger& first, const BigInteger& second) {
	int carry = 0;
		for (size_t i = 0; i < second.m_number_of_bits_int32 || carry; ++i) {
			long long tmp = first.m_bits_of_int32[i];
			tmp -= (long long)carry + (long long)(i < second.m_number_of_bits_int32 ? second.m_bits_of_int32[i] : 0);
			carry = (int)(tmp < 0);
			if (carry) tmp += two_in_32;
			first.m_bits_of_int32[i] = tmp;
		}
		first.delete_zeroes();
		return first;
}

const BigInteger operator+(BigInteger first, const BigInteger& second) {
	first += second;
	return first;
}

const BigInteger operator-(BigInteger first, const BigInteger& second) {
	first -= second;
	return first;
}

const BigInteger& BigInteger::operator++() {
	*this += 1;
	return *this;
}

const BigInteger BigInteger::operator++(int) {
	BigInteger copy (*this);
	*this += 1;
	return copy;
}

const BigInteger& BigInteger::operator--() {
	*this -= 1;
	return *this;
}
const BigInteger BigInteger::operator--(int) {
	BigInteger copy(*this);
	*this -= 1;
	return copy;
}

bool CheckString(std::string number) {
	for (size_t i = 0; i < number.size(); ++i) {
		if (number[i] < '0' || number[i]>'9') return false;
	}
	return true;
}

std::string to_string(const BigInteger& X) {
	if (X == 0) return (std::string)"0";
	BigInteger copy = X;
	if (!copy.m_sign) {
		copy.m_sign = !copy.m_sign;
	}
	std::string res="";
	do {
		unsigned int carry = 0;
		for (size_t i = 0; i < copy.m_number_of_bits_int32; ++i) {
			unsigned long long cur = copy.m_bits_of_int32[copy.m_number_of_bits_int32 - i - 1] + carry * two_in_32;
			copy.m_bits_of_int32[copy.m_number_of_bits_int32 - i - 1] = cur / 10;
			carry = cur % 10;
		}
		res += std::to_string(carry);
        if(copy.m_number_of_bits_int32>1 && copy.m_bits_of_int32[copy.m_number_of_bits_int32 - 1] == 0){
            copy.resize_down();
        }
	} while (copy != 0);
	if (!X.m_sign && !(X.m_number_of_bits_int32==1 && X.m_bits_of_int32[0]==0))
		res += '-';
	for (size_t i = 0; i < res.size() / 2; ++i) {
		std::swap(res[i], res[res.size() - i - 1]);
	}
	return res;
}

BigInteger mult_long_by_long(BigInteger& first, const BigInteger& second) {
	BigInteger res = 0;
	res.resize(first.m_number_of_bits_int32 + second.m_number_of_bits_int32);
	for (size_t i = 0; i < first.m_number_of_bits_int32; ++i) {
		long long carry = 0;
		for (size_t j = 0; j < second.m_number_of_bits_int32 || carry; ++j) {
			unsigned long long cur = (long long)res.m_bits_of_int32[i + j] + 
			(long long)first.m_bits_of_int32[i] *( j < second.m_number_of_bits_int32 ? (long long)second.m_bits_of_int32[j] : 0) 
				+ carry;
			res.m_bits_of_int32[i + j] = cur % (unsigned long long)two_in_32;
			carry = cur / two_in_32;
		}
	}
	res.delete_zeroes();
	if (!first.m_sign && second.m_sign || first.m_sign && !second.m_sign) {
		res.m_sign = false;
	}
	else {
		res.m_sign = true;
	}
	return res;
}

BigInteger& BigInteger::operator*=(const BigInteger& second){
		*this = mult_long_by_long(*this, second);
		return *this;
}

const BigInteger operator*(BigInteger first, const BigInteger& second) {
	first *= second;
	return first;
}

BigInteger BigInteger::operator!() const {
	BigInteger copy = *this;
	copy.m_sign = !copy.m_sign;
	return copy;
}

BigInteger BigInteger::operator~() const {
	BigInteger copy = *this;
	copy.m_sign = !copy.m_sign;
	--copy;
	return copy;
}

void BigInteger::additional_code() {
	if (!m_sign) {
		for (size_t i = 0; i < m_number_of_bits_int32; ++i) {
			m_bits_of_int32[i] = ~m_bits_of_int32[i];
		}
			--(*this);
	}
}

BigInteger& BigInteger::operator|=(const BigInteger& second) {
	BigInteger copy1, copy2;
	if (this->m_number_of_bits_int32 > second.m_number_of_bits_int32) {
		 copy1=*this, copy2=second;
	}
	else {
		copy2=*this, copy1=second;
	}
	copy1.additional_code();
	copy2.additional_code();
	BigInteger res(0);
	res.m_sign = (this->m_sign && second.m_sign);
	res.resize(copy1.m_number_of_bits_int32);
	for (size_t i = 0; i < copy2.m_number_of_bits_int32; ++i) {
		res.m_bits_of_int32[i] = copy1.m_bits_of_int32[i] | copy2.m_bits_of_int32[i];
	}
	if (copy2.m_sign) {
		for (size_t i = copy2.m_number_of_bits_int32; i <copy1.m_number_of_bits_int32; ++i)
			res.m_bits_of_int32[i] = copy1.m_bits_of_int32[i];
	}
	else {
		for (size_t i = copy2.m_number_of_bits_int32; i < copy1.m_number_of_bits_int32; ++i)
		res.m_bits_of_int32[i] = 1;
	}
	res.delete_zeroes();
	res.additional_code();
	*this = res;
	return *this;
}

const BigInteger operator|(BigInteger first, const BigInteger& second) {
	first |= second;
	return first;
}

BigInteger& BigInteger::operator&=(const BigInteger& second) {
	BigInteger copy1, copy2;
	if (this->m_number_of_bits_int32 > second.m_number_of_bits_int32) {
		copy1 = *this;
		copy2 = second;
	}
	else {
		copy2 = *this;
		copy1 = second;
	}
	copy1.additional_code();
	copy2.additional_code();
	BigInteger res(0);
	res.m_sign = (this->m_sign || second.m_sign);
	res.resize(copy1.m_number_of_bits_int32);
	for (size_t i = 0; i < copy2.m_number_of_bits_int32; ++i) {
		res.m_bits_of_int32[i] = copy1.m_bits_of_int32[i] & copy2.m_bits_of_int32[i];
	}
	if (copy2.m_sign) {
		for (size_t i = copy2.m_number_of_bits_int32; i < copy1.m_number_of_bits_int32; ++i)
			res.m_bits_of_int32[i] = copy1.m_bits_of_int32[i];
	}
	else {
		for (size_t i = copy2.m_number_of_bits_int32; i < copy1.m_number_of_bits_int32; ++i)
			res.m_bits_of_int32[i] = 1;
	}
	res.delete_zeroes();
	res.additional_code();
	*this = res;
	return *this;
}

const BigInteger operator&(BigInteger first, const BigInteger& second) {
	first &= second;
	return first;
}

BigInteger& BigInteger::operator^=(const BigInteger& second) {
	BigInteger copy1, copy2;
	if (this->m_number_of_bits_int32 > second.m_number_of_bits_int32) {
		copy1 = *this;
		copy2 = second;
	}
	else {
		copy2 = *this;
		copy1 = second;
	}
	copy1.additional_code();
	copy2.additional_code();
	BigInteger res(0);
	res.m_sign = !(this->m_sign ^ second.m_sign);
	res.resize(copy1.m_number_of_bits_int32);
	for (size_t i = 0; i < copy2.m_number_of_bits_int32; ++i) {
		res.m_bits_of_int32[i] = copy1.m_bits_of_int32[i] ^ copy2.m_bits_of_int32[i];
	}
	if (copy2.m_sign) {
		for (size_t i = copy2.m_number_of_bits_int32; i < copy1.m_number_of_bits_int32; ++i)
			res.m_bits_of_int32[i] = copy1.m_bits_of_int32[i];
	}
	else {
		for (size_t i = copy2.m_number_of_bits_int32; i < copy1.m_number_of_bits_int32; ++i)
			res.m_bits_of_int32[i] = 1;
	}
	res.delete_zeroes();
	res.additional_code();
	*this = res;
	return *this;
}

const BigInteger operator^(BigInteger first, const BigInteger& second) {
	first ^= second;
	return first;
}

BigInteger module(BigInteger object) {
	if (object < 0) {
		return -object;
	}
	return object;
}

BigInteger& BigInteger::operator/=(const BigInteger& second) {
	if (module(*this) < module(second) && second!=0 || *this==0) {
		*this = 0;
		return *this;
	}
	if (second.m_number_of_bits_int32 == 1) {
		if (second.m_bits_of_int32[0] != 0) {
			*this = div_long_by_short(*this, second);
			this->m_sign = (this->m_sign && second.m_sign || !this->m_sign && !second.m_sign);
			return *this;
		}
		else {
			throw std::invalid_argument(to_string(second));
		}
	}
	BigInteger two(2), copy1(*this), copy2(second);
	copy1.m_sign = true;
	copy2.m_sign = true;
	BigInteger right, left, middle = div_long_by_short(copy1, two), tmp, prev_middle;
	while(middle!=prev_middle){
		if (middle * copy2 > copy1) {
			tmp = middle;
			prev_middle = middle;
			middle = div_long_by_short(middle + left, two);
			right = tmp;
		}
		else if (middle * copy2 < copy1) {
			tmp = middle;
			prev_middle = middle;
			middle = div_long_by_short(middle + right, two);
			left = tmp;
		}
		else {
			break;
		}
	}
	bool copy_sign = this->m_sign;
	*this = middle;
	this->m_sign = (copy_sign && second.m_sign || !copy_sign && !second.m_sign);
	return *this;
}

const BigInteger operator/(BigInteger first, const BigInteger& second) {
	first /= second;
	return first;
}

BigInteger div_long_by_short(BigInteger first, const BigInteger second) {
	unsigned int carry = 0;
	for (int i = (int)first.m_number_of_bits_int32 - 1; i >= 0; --i) {
		unsigned long long cur = first.m_bits_of_int32[i] + carry * 1ll *two_in_32;
		first.m_bits_of_int32[i] = (unsigned int)(cur / second.m_bits_of_int32[0]);
		carry = (unsigned int)(cur % second.m_bits_of_int32[0]);
	}
	first.delete_zeroes();
	return first;
}

BigInteger& BigInteger::operator%=(const BigInteger& second) {
	*this = *this - (*this / second) * second;
	return *this;
}

BigInteger operator%(BigInteger first, const BigInteger& second) {
	first %= second;
	return first;
}

BigInteger& BigInteger::operator>>=(int step) {
	if (step < 0)
		throw std::invalid_argument("");
	size_t difference_between_cur_and_res = step / 32;
	if (difference_between_cur_and_res > this->m_number_of_bits_int32) {
		*this = 0;
		return *this;
	}
	size_t number_of_passing_bits = step % 32;
	unsigned int remain = 1 << (number_of_passing_bits);
	unsigned int remain_next = 0;
	unsigned int remain_this;
	this->m_number_of_bits_int32 -= difference_between_cur_and_res;
	for (int i = this->m_number_of_bits_int32-1; i >-1; --i) {
		remain_this = remain_next;
		remain_next = this->m_bits_of_int32[i] % remain;
		this->m_bits_of_int32[i] >>= number_of_passing_bits;
		this->m_bits_of_int32[i] |= (remain_this << (32 - number_of_passing_bits));
	}
	unsigned int* tmp = new unsigned int[this->m_number_of_bits_int32];
	if (difference_between_cur_and_res) {
		for (int i = this->m_number_of_bits_int32-1; i >-1; --i)
			tmp[i] = this->m_bits_of_int32[i];
		delete[] m_bits_of_int32;
		this->m_bits_of_int32 = tmp;
		tmp = nullptr;
		delete[] tmp;
	}
	this->delete_zeroes();
	if (!m_sign) {
		m_sign = true;
		*this = ~(*this);
	}
	return *this;
}

const BigInteger operator>>(BigInteger first, int step) {
	first >>= step;
	return first;
}

BigInteger& BigInteger::operator<<=(int step) {
	if (step < 0)
		throw std::invalid_argument("");
	size_t new_bits_after_shift = step / 32;
	size_t rest = step % 32;
	if (new_bits_after_shift)
		m_number_of_bits_int32 += new_bits_after_shift;
	if (rest>0) 
		++m_number_of_bits_int32;
	unsigned int* tmp = new unsigned int[m_number_of_bits_int32];
	for (size_t i = 0; i < new_bits_after_shift; ++i)
		tmp[i] = 0;
	for (size_t i = new_bits_after_shift; i < m_number_of_bits_int32; ++i)
		tmp[i] = m_bits_of_int32[i - new_bits_after_shift];
	if (rest)
		tmp[m_number_of_bits_int32-1] = 0;
	unsigned int x = two_in_32 - 1;
	unsigned int remain_next = 0;
	unsigned int remain_this;
	if (rest) {
		for (size_t i = 0; i < m_number_of_bits_int32; ++i) {
			remain_this = remain_next;
			remain_next = tmp[i] | (x >> rest);
			tmp[i] <<= rest;
			tmp[i] |= (remain_this >> (32 - rest));
		}
	}
	delete[] m_bits_of_int32;
	m_bits_of_int32 = tmp;
	tmp = nullptr;
	delete[] tmp;
	this->delete_zeroes();
	return *this;
}

const BigInteger operator<<(BigInteger first, int step) {
	first <<= step;
	return first;
}


