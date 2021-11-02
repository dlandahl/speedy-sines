
/*

Status:

Scalar with lookup table:
   ~0.89s, sounds good

Scalar with -ffast-math sin:
   ~1.70s, sounds good

Scalar with nick_sine:
   ~0.80 with mild spectral artifacting

Vector with SSE_MATHFUN:
   ~1.12s, sounds good

Vector with AVX_MATHFUN:
   ~0.26s sounds good

Vector with AVX Gather lookup table:
   ~0.16s sounds good

*/

#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>

#include <stdint.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#define WIN32_LEAN_AND_MEAN
#include <windows.h>

#include <immintrin.h>

#if defined(__AVX__)
  #define AVX_ENABLED 1
#else
  #define AVX_ENABLED 0
#endif

#if defined(__AVX2__)
  #define AVX2_ENABLED 1
#else
  #define AVX2_ENABLED 0
#endif

#if defined(__FAST_MATH__)
  #define FASTMATH_ENABLED 1
#else
  #define FASTMATH_ENABLED 0
#endif

#define SSE_MATHFUN_WITH_CODE
#include "sse_mathfun/sse_mathfun.h"

#if AVX_ENABLED
  #include "avx_mathfun/avx_mathfun.h"
#endif

typedef uint32_t u32;
typedef uint16_t u16;
typedef uint8_t  u8;
typedef float    f32;
typedef double   f64;

#define WAVE_ID(text) \
(u32) (text[0] | (text[1] << 8) | (text[2] << 16) | (text[3] << 24))

#define SAMPLE_RATE 48000
#define tau 6.28318530717958647693
#define pi  3.14159265358979323846

#define TABLE_SIZE 512
f32 sin_table[TABLE_SIZE];

void generate_sine_table() {
    for (int n = 0; n < TABLE_SIZE; n++) {
        sin_table[n] = sin(n * tau / TABLE_SIZE);
    }
}

f32 read(f32 phase) {
    phase *= TABLE_SIZE;
    int a = phase;
    int b = phase + 1;
    f32 t = phase - a;

    return sin_table[a] * (1 - t) + sin_table[b] * t;
}

#define FORMAT_PCM        0x0001;
#define FORMAT_IEEE_FLOAT 0x0003;
#define FORMAT_ALAW       0x0006;
#define FORMAT_MULAW      0x0007;
#define FORMAT_EXTENSIBLE 0xFFFE;

struct Riff_Header {
    u32 header_id;
    u32 header_size;
    u32 wave_id;
    u32 fmt_id;
    u32 fmt_size;
    u16 fmt_tag;
    u16 number_of_channels;
    u32 sample_rate;
    u32 bytes_per_second;
    u16 block_align;
    u16 bit_depth;
    u32 data_id;
    u32 data_size;
};

struct Riff_Header create_header(int number_of_samples) {
    int bytes_per_sample = sizeof(f32);
    int number_of_channels = 1;

    struct Riff_Header header;

    header.header_id          = WAVE_ID("RIFF");
    header.header_size        = 36 + number_of_samples * bytes_per_sample;
    header.wave_id            = WAVE_ID("WAVE");
    header.fmt_id             = WAVE_ID("fmt ");
    header.fmt_size           = 16;
    header.fmt_tag            = FORMAT_IEEE_FLOAT;
    header.number_of_channels = number_of_channels;
    header.sample_rate        = SAMPLE_RATE;
    header.bytes_per_second   = SAMPLE_RATE * bytes_per_sample * number_of_channels;
    header.block_align        =               bytes_per_sample * number_of_channels;
    header.bit_depth          = bytes_per_sample * 8;
    header.data_id            = WAVE_ID("data");
    header.data_size          = number_of_samples * bytes_per_sample;

    return header;
}

void write_wav_file(f32* data, int count, const char* filename) {
    int header_size = sizeof(struct Riff_Header);
    assert(header_size == 44);
    u8* buffer = malloc(sizeof(f32) * count + header_size);

    struct Riff_Header header = create_header(count);
    *(struct Riff_Header*) buffer = header;
    memcpy(buffer + header_size, data, sizeof(f32) * count);

    FILE* file = fopen(filename, "wb");
    int result = fwrite(buffer, 1, sizeof(f32) * count + header_size, file);

    fclose(file);
    free(buffer);
}

// 1seconds = 48000samples
// 1cycles  = 48000samples / hz


void table_scalar_additive_synthesis(f32* data, int count) {
    memset(data, 0, count * sizeof(f32));
    f64 cycles_per_sample = 50.0 / SAMPLE_RATE;

    for (int p = 1; p <= 10; p++) {
        f64 phase = 0;
        for (int n = 0; n < count; n++) {
            data[n] += read(phase) / p;
            phase += (f64) p * cycles_per_sample;
            if (phase >= 1) phase -= 1;
        }
    }
}

// Sine function by "Nick".
// http://web.archive.org/web/20141220225551/http://forum.devmaster.net/t/fast-and-accurate-sine-cosine/9648
float nick_sine(f32 x) {
    f32 b =  4 /  pi;
    f32 c = -4 / (pi*pi);

    f32 y = b * x + c * x * fabsf(x);

    f32 p = 0.225;
    return p * (y * fabsf(y) - y) + y;
}

void nick_scalar_additive_synthesis(f32* data, int count) {
    memset(data, 0, count * sizeof(f32));
    f64 cycles_per_sample = tau * 50.0 / SAMPLE_RATE;

    for (int p = 1; p <= 10; p++) {
        f64 phase = 0;
        for (int n = 0; n < count; n++) {
            data[n] += nick_sine(phase) / p;
            phase += (f64) p * cycles_per_sample;
            phase -= tau * (phase > pi);
        }
    }
}

void fastmath_scalar_additive_synthesis(f32* data, int count) {
    memset(data, 0, count * sizeof(f32));
    f64 cycles_per_sample = tau * 50.0 / SAMPLE_RATE;

    for (int p = 1; p <= 10; p++) {
        f64 phase = 0;
        for (int n = 0; n < count; n++) {
            data[n] += sin(phase) / p;
            phase += (f64) p * cycles_per_sample;
            if (phase > pi) phase -= tau;
        }
    }
}

void vector_additive_synthesis(f32* data, int count) {
    memset(data, 0, count * sizeof(f32));

    assert(count % 4 == 0);
    f64 mul = tau * 50.0 / SAMPLE_RATE;

    __m128 indices = _mm_setr_ps(0, 1, 2, 3);
    __m128 vtau    = _mm_set1_ps(tau);

    for (int p = 1; p <= 10; p++) {
        f64 phase = 0;
        f64 step_size = p * mul;

        __m128 amplitude = _mm_set1_ps(1.0 / p);
        __m128 offset = _mm_mul_ps(_mm_set1_ps(step_size), indices);

        for (int n = 0; n < count; n += 4) {
            __m128 vdata  = _mm_load_ps(data + n);
            __m128 vphase = _mm_add_ps(_mm_set1_ps(phase), offset);

            __m128 mask = _mm_cmpge_ps(vphase, vtau);
            vphase      = _mm_sub_ps  (vphase, _mm_and_ps(vtau, mask));

            __m128 sine =     sin_ps(vphase);
            sine        = _mm_mul_ps(sine, amplitude);
            vdata       = _mm_add_ps(sine, vdata);
            _mm_store_ps(data + n, vdata);

            phase += step_size * 4;
            phase -= tau * (phase > tau);
        }
    }
}

#if AVX_ENABLED
void vector256_additive_synthesis(f32* data, int count) {
    memset(data, 0, count * sizeof(f32));

    assert(count % 8 == 0);
    f64 mul = tau * 50.0 / SAMPLE_RATE;

    __m256 indices = _mm256_setr_ps(0, 1, 2, 3, 4, 5, 6, 7);
    __m256 max     = _mm256_set1_ps(tau);

    for (int p = 1; p <= 10; p++) {
        f64 phase = 0;
        f64 step_size = p * mul;

        __m256 amplitude = _mm256_set1_ps(1.0 / p);
        __m256 offset = _mm256_mul_ps(_mm256_set1_ps(step_size), indices);

        for (int n = 0; n < count; n += 8) {
            __m256 vdata  = _mm256_load_ps(data + n);
            __m256 vphase = _mm256_add_ps(_mm256_set1_ps(phase), offset);

            __m256 mask = _mm256_cmp_ps(vphase, max, _CMP_NLT_UQ);
            vphase      = _mm256_sub_ps(vphase, _mm256_and_ps(max, mask));

            __m256 sine =     sin256_ps(vphase);
            sine        = _mm256_mul_ps(sine, amplitude);
            vdata       = _mm256_add_ps(sine, vdata);
            _mm256_store_ps(data + n, vdata);

            phase += step_size * 8;
            phase -= tau * (phase > tau);
        }
    }
}
#endif

#if AVX2_ENABLED
__m256 read_vector(__m256 phase) {
    phase = _mm256_mul_ps(phase, _mm256_set1_ps(TABLE_SIZE));
    __m256i a = _mm256_cvtps_epi32(_mm256_floor_ps(phase));
    __m256i b = _mm256_cvtps_epi32(_mm256_floor_ps(_mm256_add_ps(phase, _mm256_set1_ps(1))));

    __m256 t     = _mm256_sub_ps(phase, _mm256_floor_ps(phase));
    __m256 lower = _mm256_i32gather_ps(sin_table, a, 4);
    __m256 upper = _mm256_i32gather_ps(sin_table, b, 4);

    upper = _mm256_mul_ps(upper, t);
    t     = _mm256_sub_ps(_mm256_set1_ps(1), t);
    lower = _mm256_mul_ps(lower, t);
    return _mm256_add_ps(lower, upper);
}

void vector256_table_additive_synthesis(f32* data, int count) {
    memset(data, 0, count * sizeof(f32));

    assert(count % 8 == 0);
    f64 mul = 440.0 / SAMPLE_RATE;

    __m256 indices = _mm256_setr_ps(0, 1, 2, 3, 4, 5, 6, 7);
    __m256 max     = _mm256_set1_ps(1);

    for (int p = 1; p <= 10; p++) {
        f64 phase = 0;
        f64 step_size = p * mul;

        __m256 amplitude = _mm256_set1_ps(1.0 / p);
        __m256 offset = _mm256_mul_ps(_mm256_set1_ps(step_size), indices);

        for (int n = 0; n < count; n += 8) {
            __m256 vdata  = _mm256_load_ps(data + n);
            __m256 vphase = _mm256_add_ps(_mm256_set1_ps(phase), offset);

            __m256 mask = _mm256_cmp_ps(vphase, max, _CMP_NLT_UQ);
            vphase      = _mm256_sub_ps(vphase, _mm256_and_ps(max, mask));

            __m256 sine =  read_vector(vphase);
            vdata       = _mm256_add_ps(sine, vdata);
            _mm256_store_ps(data + n, vdata);

            phase += step_size * 8;
            while (phase > 1) {
                phase -= 1;
            }
      }
}
#endif

#define REPEAT 1
#define PROFILE(stmt, id) {                                                      \
    LARGE_INTEGER before;                                                        \
    LARGE_INTEGER after;                                                         \
    f64 elapsed = 0;                                                             \
    for (int n = 0; n < REPEAT; n++) {                                           \
        QueryPerformanceCounter(&before);                                        \
        stmt;                                                                    \
        QueryPerformanceCounter(&after);                                         \
        elapsed += (f64) (after.QuadPart - before.QuadPart) * tick_period;       \
    }                                                                            \
    printf("| %lfs | %s\n", elapsed / REPEAT, #id);                              \
}

int main() {
    if (!FASTMATH_ENABLED) {
        printf("Error! Compile with -ffast-math.\n");
    }
    generate_sine_table();

#if AVX2_ENABLED
    printf("Running with AVX2 instructions enabled.\n");
#elif AVX_ENABLED
    printf("Running with AVX instructions enabled.\n");
#else
    printf("AVX instructions are disabled.\n");
#endif

    LARGE_INTEGER performance_frequency;
    QueryPerformanceFrequency(&performance_frequency);
    f64 tick_period = 1.0 / (f64) performance_frequency.QuadPart;

    int count = SAMPLE_RATE * 60 * 10;
    f32* data = malloc(sizeof(f32) * count);

    PROFILE(table_scalar_additive_synthesis(data, count), SCALAR_TABLE);
    write_wav_file(data, count, "scalar_table.wav");
    
    PROFILE(fastmath_scalar_additive_synthesis(data, count), SCALAR_FASTMATH);
    write_wav_file(data, count, "scalar_fastmath.wav");
    
    PROFILE(nick_scalar_additive_synthesis(data, count), SCALAR_CUSTOM_FAST_SINE);
    write_wav_file(data, count, "scalar_custom.wav");
    
    PROFILE(vector_additive_synthesis(data, count), SSE_MATHFUN);
    write_wav_file(data, count, "sse_mathfun.wav");

#if AVX_ENABLED
    PROFILE(vector256_additive_synthesis(data, count), AVX_MATHFUN);
    write_wav_file(data, count, "avx_mathfun.wav");
#endif

#if AVX2_ENABLED
    PROFILE(vector256_table_additive_synthesis(data, count), AVX2_TABLE);
    write_wav_file(data, count, "avx2_table.wav");
#endif

    free(data);
    return 0;
}
