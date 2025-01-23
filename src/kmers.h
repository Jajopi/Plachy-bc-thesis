#pragma once

#include <string>
#include <iostream>
#include <cstdint>

#include "uint256_t/uint256_t.h"

#include "ac/kmers_ac.h"

typedef __uint128_t kmer128_t;
typedef uint64_t kmer64_t;
typedef uint256_t kmer256_t;

static const uint8_t nucleotideToInt[] = {
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
	};

/// Compute the prefix of size d of the given k-mer.
template <typename kmer_t>
kmer_t BitPrefix(kmer_t kMer, int k, int d) {
    return kMer >> ((k - d) << kmer_t(1));
}

/// Compute the suffix of size d of the given k-mer.
template <typename kmer_t>
kmer_t BitSuffix(kmer_t kMer, int d) {
    return kMer & ((kmer_t(1) << (d << kmer_t(1))) - kmer_t(1));
}

/// Checkered mask. cmask<uint16_t, 1> is every other bit on
/// (0x55). cmask<uint16_t,2> is two bits one, two bits off (0x33). Etc.
/// Copyright: Jellyfish GPL-3.0
template<typename U, int len, int l = sizeof(U) * 8 / (2 * len)>
struct cmask {
    static const U v =
            (cmask<U, len, l - 1>::v << (2 * len)) | (((U)1 << len) - 1);
};
template<typename U, int len>
struct cmask<U, len, 0> {
    static const U v = 0;
};

/// Compute the reverse complement of a word.
/// Copyright: Jellyfish GPL-3.0
inline kmer64_t word_reverse_complement(kmer64_t w) {
    typedef kmer64_t U;
    w = ((w >> 2)  & cmask<U, 2 >::v) | ((w & cmask<U, 2 >::v) << 2);
    w = ((w >> 4)  & cmask<U, 4 >::v) | ((w & cmask<U, 4 >::v) << 4);
    w = ((w >> 8)  & cmask<U, 8 >::v) | ((w & cmask<U, 8 >::v) << 8);
    w = ((w >> 16) & cmask<U, 16>::v) | ((w & cmask<U, 16>::v) << 16);
    w = ( w >> 32                   ) | ( w                    << 32);
    return ((U)-1) - w;
}

/// Compute the reverse complement of a word.
/// Copyright: Jellyfish GPL-3.0
inline kmer128_t word_reverse_complement(kmer128_t w) {
    typedef kmer128_t U;
    w = ((w >> 2)  & cmask<U, 2 >::v) | ((w & cmask<U, 2 >::v) << 2);
    w = ((w >> 4)  & cmask<U, 4 >::v) | ((w & cmask<U, 4 >::v) << 4);
    w = ((w >> 8)  & cmask<U, 8 >::v) | ((w & cmask<U, 8 >::v) << 8);
    w = ((w >> 16) & cmask<U, 16>::v) | ((w & cmask<U, 16>::v) << 16);
    w = ((w >> 32) & cmask<U, 32>::v) | ((w & cmask<U, 32>::v) << 32);
    w = ( w >> 64                   ) | ( w                    << 64);
    return ((U)-1) - w;
}

/// Compute the reverse complement of a word.
inline kmer256_t word_reverse_complement(kmer256_t w) {
    kmer128_t low = word_reverse_complement(w.lower());
    kmer128_t high = word_reverse_complement(w.upper());
    return kmer256_t(low, high);
}

/// Compute the reverse complement of the given k-mer.
template <typename kmer_t>
kmer_t ReverseComplement(kmer_t kMer, int k) {
    return (((kmer_t)word_reverse_complement(kMer)) >> ((sizeof(kMer)<<3) - (k << 1))) & ((kmer_t(1) << (k << 1)) - kmer_t(1));
}

/// Return the number 0 - 3 correcponding to index-th nucleotide from the encoded k-mer.
template <typename kmer_t>
inline uint8_t NucleotideIndexAtIndex(kmer_t encoded, int k, int index) {
    return (uint8_t)((encoded >> ((k - index - kmer_t(1)) << kmer_t(1))) & kmer_t(3));
}
template <typename kmer_t, typename other>
inline uint8_t NucleotideIndexAtIndex(std::pair<kmer_t, other> p, int k, int index) {
    return (uint8_t)((p.first >> ((k - index - kmer_t(1)) << kmer_t(1))) & kmer_t(3));
}

const char letters[4] {'A', 'C', 'G', 'T'};
/// Return the index-th nucleotide from the encoded k-mer.
template <typename kmer_t>
inline char NucleotideAtIndex(kmer_t encoded, int k, int index) {
    return letters[NucleotideIndexAtIndex(encoded, k, index)];
}
template <typename kmer_t, typename other>
inline char NucleotideAtIndex(std::pair<kmer_t, other> p, int k, int index) {
    return letters[NucleotideIndexAtIndex(p, k, index)];
}

/// Convert the encoded KMer representation to string.
template <typename kmer_t>
std::string NumberToKMer(kmer_t encoded, int length) {
    std::string ret(length, 'N');
    for (int i = 0; i < length; ++i) {
        // The last two bits correspond to one nucleotide.
        ret[length - i -1] = letters[(uint64_t)(encoded & 3)];
        // Move to the next letter.
        encoded >>= 2;
    }
    return ret;
}

// Happily reimplementing the wheel
char to_lower(char c){
    return c | 32;
}
char to_upper(char c){
    return c & 95;
}

template <typename kmer_t, typename size_t_max>
void print_kmer_masked(kmer_t kmer, size_t_max k, std::ostream& os, size_t_max prefix = 0){
    if (prefix == 0 || prefix > k) prefix = k;

    os << to_upper(NucleotideAtIndex(kmer, k, 0));
    for (size_t_max c = 1; c < prefix; ++c){
        os << to_lower(NucleotideAtIndex(kmer, k, c));
    }
}

template <typename kmer_t, typename size_t_max>
void print_kmer(kmer_t kmer, size_t_max k, std::ostream& os, size_t_max prefix = 0){
    if (prefix == 0 || prefix > k) prefix = k;

    for (size_t_max c = 0; c < prefix; ++c){
        os << to_upper(NucleotideAtIndex(kmer, k, c));
    }
}

template<typename kmer_t, typename size_t_max>
size_t_max compute_max_overlap(kmer_t kmer1, kmer_t kmer2, size_t_max k){
    for (size_t_max ov = k; ov > 0; --ov){
        bool found = true;
        for (size_t_max i = 0; i < ov; ++i){
            if (NucleotideAtIndex(kmer2, k, i) != NucleotideAtIndex(kmer1, k, k - ov + i)){
                found = false;
                break;
            }
        }
        if (found) return ov;
    }
    return 0;
}

template <typename kmer_t, typename size_t_max>
size_t_max decode_and_print_indexes(const std::vector<kmer_t>& kMers, const std::vector<size_t_max>& indexes, std::ostream& os, size_t_max k,
    bool encode_mask = true, bool count_not_print = false){
    size_t_max total_length = 0;
    size_t_max run_count = 1;
    
    kmer_t actual_kmer = kMers[indexes[0]], new_kmer = 0;
    for (size_t_max i = 1; i < indexes.size(); ++i){
        new_kmer = kMers[indexes[i]];
        
        size_t_max ov = compute_max_overlap(actual_kmer, new_kmer, k);
        //print_kmer(actual_kmer, k, os, size_t_max(k)); std::cout << ' '; print_kmer(new_kmer, k, os, size_t_max(k)); std::cout << ' ' << ov << std::endl;
        if (encode_mask && !count_not_print) print_kmer_masked(actual_kmer, k, os, size_t_max(k - ov));
        else if(!count_not_print) print_kmer(actual_kmer, k, os, size_t_max(k - ov));
        
        total_length += k - ov;
        if (ov < k - 1) ++run_count;

        actual_kmer = new_kmer;
    }
    if (encode_mask && !count_not_print) print_kmer_masked(new_kmer, k, os);
    else if(!count_not_print) print_kmer(new_kmer, k, os);
    
    total_length += k;

    os << std::endl;
    std::cerr << "Total length: " << total_length << std::endl;
    // std::cerr << "Run count: " << run_count << std::endl;

    return total_length;
}
