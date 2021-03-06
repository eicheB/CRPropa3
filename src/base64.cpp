/// base64 encodig and decoding.
///
/// Based on the implementations by
/// Jouni Malinen <j@w1.fi> and contributors from wpa_supplicant and hostapd in
/// http://web.mit.edu/freebsd/head/contrib/wpa/ and
/// http://web.mit.edu/freebsd/head/contrib/wpa/src/utils/base64.c and
/// http://web.mit.edu/freebsd/head/contrib/wpa/src/utils/base64.h
///
/// Published under a 3-clause BSD license
///
#include <string>
#include <cstring>
#include <stdexcept>

#include "crpropa/base64.h"

namespace crpropa 
{

//Alphabet used
static const unsigned char encode_alphabet[65] =
"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
static int decode_alphabet[256] = {-1};

/// Encodes data
std::string Base64::encode(const unsigned char *src, size_t len)
{
		size_t olen = 4*((len + 2) / 3); /* 3-byte blocks to 4-byte */

		if (olen < len)
			throw std::runtime_error("Integer overflow in Base64::encoding, data to large!");

		std::string outStr;
		outStr.resize(olen);

		unsigned char *out = (unsigned char*) outStr.c_str();
		unsigned char *pos = out;
		const unsigned char *end, *in;

		end = src + len;
		in = src;
		while (end - in >= 3) {
				*pos++ = encode_alphabet[in[0] >> 2];
				*pos++ = encode_alphabet[((in[0] & 0x03) << 4) | (in[1] >> 4)];
				*pos++ = encode_alphabet[((in[1] & 0x0f) << 2) | (in[2] >> 6)];
				*pos++ = encode_alphabet[in[2] & 0x3f];
				in += 3;
		}

		if (end - in) {
				*pos++ = encode_alphabet[in[0] >> 2];
				if (end - in == 1) {
						*pos++ = encode_alphabet[(in[0] & 0x03) << 4];
						*pos++ = '=';
				}
				else {
						*pos++ = encode_alphabet[((in[0] & 0x03) << 4) |
								(in[1] >> 4)];
						*pos++ = encode_alphabet[(in[1] & 0x0f) << 2];
				}
				*pos++ = '=';
		}

		return outStr;
}


std::string Base64::decode(const std::string &data)
{
const unsigned char *src = (unsigned char*) data.c_str();
size_t len = data.size();

if (decode_alphabet[0] == -1)
{ // build decode alphabet
	std::memset(decode_alphabet, 0x80, 256);
	for (size_t i = 0; i < sizeof(encode_alphabet) - 1; i++)
		decode_alphabet[encode_alphabet[i]] = (unsigned char) i;
	decode_alphabet['='] = 0;
}

size_t olen = 0;
for (size_t i = 0; i < len; i++) {
	if (decode_alphabet[src[i]] != 0x80)
		olen++;
}

if (olen == 0 || olen % 4)
			throw std::runtime_error("Base64 decode, invalid input size");

olen = olen / 4 * 3;
std::string str;
str.resize(olen);

unsigned char *out = (unsigned char*) str.c_str();
unsigned char *pos = out;

size_t count = 0;
int pad = 0;
unsigned char block[4];
for (size_t i = 0; i < len; i++) {
	unsigned char tmp = decode_alphabet[src[i]];
	if (tmp == 0x80)
		continue;

	if (src[i] == '=')
		pad++;
	block[count] = tmp;
	count++;
	if (count == 4) {
		*pos++ = (block[0] << 2) | (block[1] >> 4);
		*pos++ = (block[1] << 4) | (block[2] >> 2);
		*pos++ = (block[2] << 6) | block[3];
		count = 0;
		if (pad) {
			if (pad == 1)
				pos--;
			else if (pad == 2)
				pos -= 2;
			else {
				throw std::runtime_error("Base64 decode, invalid padding");
			}
			break;
		}
	}
}
return str;
}
};
