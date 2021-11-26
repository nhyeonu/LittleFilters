#ifndef LITTLEFILTERS_H_
#define LITTLEFILTERS_H_

#ifdef __cplusplus
extern "C" {
#endif

#define LF_FLOATING_POINT float

typedef struct lf_filter lf_filter;
lf_filter* lf_filter_init(const LF_FLOATING_POINT *impulse_response, int impulse_response_length, int processing_chunk_length);
void lf_filter_uninit(lf_filter* filter);

void lf_filter_replace_impulse_response(lf_filter *filter, const LF_FLOATING_POINT *impulse_response);
void lf_filter_process(lf_filter* filter, const LF_FLOATING_POINT* input, LF_FLOATING_POINT* output);

typedef enum { LF_BLACKMAN_WINDOW, LF_HAMMING_WINDOW } lf_window_func_t;
void lf_windowed_sinc_lowpass(LF_FLOATING_POINT *output, int N, lf_window_func_t window_func, LF_FLOATING_POINT cutoff_freq);
void lf_windowed_sinc_highpass(LF_FLOATING_POINT *output, int N, lf_window_func_t window_func, LF_FLOATING_POINT cutoff_freq);
void lf_windowed_sinc_bandpass(LF_FLOATING_POINT *output, int N, lf_window_func_t window_func, LF_FLOATING_POINT cutoff_freq, LF_FLOATING_POINT band_width);
void lf_windowed_sinc_bandreject(LF_FLOATING_POINT *output, int N, lf_window_func_t window_func, LF_FLOATING_POINT cutoff_freq, LF_FLOATING_POINT band_width);

void lf_moving_average(LF_FLOATING_POINT *output, int N);

void lf_fft(const LF_FLOATING_POINT *in_real, const LF_FLOATING_POINT *in_imag, LF_FLOATING_POINT *out_real, LF_FLOATING_POINT *out_imag, int N);
void lf_ifft(const LF_FLOATING_POINT *in_real, const LF_FLOATING_POINT *in_imag, LF_FLOATING_POINT *out_real, LF_FLOATING_POINT *out_imag, int N);

#ifdef __cplusplus
}
#endif

#endif
