#include "littlefilters.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>

#define LF_E 2.718281828459045

void lf_complex_multiply(LF_FLOATING_POINT a_real, LF_FLOATING_POINT a_imag,
                         LF_FLOATING_POINT b_real, LF_FLOATING_POINT b_imag,
                         LF_FLOATING_POINT *out_real, LF_FLOATING_POINT *out_imag)
{
    *out_real = a_real * b_real - a_imag * b_imag;
    *out_imag = a_real * b_imag + a_imag * b_real;
}

void lf_complex_exp(LF_FLOATING_POINT x_real, LF_FLOATING_POINT x_imag, 
                    LF_FLOATING_POINT *out_real, LF_FLOATING_POINT *out_imag)
{
    lf_complex_multiply(pow(LF_E, x_real), 0, cos(x_imag), sin(x_imag), out_real, out_imag);
}

void lf_fft_recersion(const LF_FLOATING_POINT *in_real, const LF_FLOATING_POINT *in_imag, 
                      LF_FLOATING_POINT *out_real, LF_FLOATING_POINT *out_imag, 
                      int N, int increment)
{
    if(N == 1)
    {
        *out_real = *in_real;
        *out_imag = *in_imag;
    }
    else
    {
        lf_fft_recersion(in_real, in_imag, out_real, out_imag, N / 2, increment * 2);
        lf_fft_recersion(in_real + increment, in_imag + increment, out_real + N / 2, out_imag + N / 2, N / 2, increment * 2);

        for(int k = 0; k < N/2; k++)
        {
            LF_FLOATING_POINT even_k_real = out_real[k];
            LF_FLOATING_POINT even_k_imag = out_imag[k];

            LF_FLOATING_POINT exp_real, exp_imag;
            lf_complex_exp(0, -2 * M_PI / N * k, &exp_real, &exp_imag);

            LF_FLOATING_POINT multiple_real, multiple_imag;
            lf_complex_multiply(exp_real, exp_imag, out_real[k + N / 2], out_imag[k + N / 2], &multiple_real, &multiple_imag);

            out_real[k] = even_k_real + multiple_real;
            out_imag[k] = even_k_imag + multiple_imag;
            out_real[k + N / 2] = even_k_real - multiple_real;
            out_imag[k + N / 2] = even_k_imag - multiple_imag;
        }
    }
}

void lf_fft(const LF_FLOATING_POINT *in_real, const LF_FLOATING_POINT *in_imag, LF_FLOATING_POINT *out_real, LF_FLOATING_POINT *out_imag, int N)
{
    lf_fft_recersion(in_real, in_imag, out_real, out_imag, N, 1);
}

void lf_ifft(const LF_FLOATING_POINT *in_real, const LF_FLOATING_POINT *in_imag, LF_FLOATING_POINT *out_real, LF_FLOATING_POINT *out_imag, int N)
{
    lf_fft(in_imag, in_real, out_imag, out_real, N);

    for(int i = 0; i < N; i++)
    {
        out_real[i] /= N;
        out_imag[i] /= N;
    }
}

LF_FLOATING_POINT lf_blackman_window(LF_FLOATING_POINT x, int N)
{
    return 0.42 - 0.5 * cos(2 * M_PI * x / (N - 1)) + 0.08 * cos(4 * M_PI * x / (N - 1));
}

LF_FLOATING_POINT lf_hamming_window(LF_FLOATING_POINT x, int N)
{
    return 0.54 - 0.46 * cos(2 * M_PI * x / (N - 1));
}

LF_FLOATING_POINT lf_normalized_sinc(LF_FLOATING_POINT x)
{
    if(x == 0)
        return 1;
    else
        return sin(M_PI * x) / (M_PI * x);
}

LF_FLOATING_POINT lf_windowed_sinc(LF_FLOATING_POINT x, LF_FLOATING_POINT cutoff_frequency, lf_window_func_t window_func, int filter_size)
{
    switch(window_func)
    {
        case LF_BLACKMAN_WINDOW:
            return 2 * cutoff_frequency * lf_normalized_sinc(2 * cutoff_frequency * (x - (filter_size - 1) / 2.0)) * lf_blackman_window(x, filter_size - 1);
        case LF_HAMMING_WINDOW:
            return 2 * cutoff_frequency * lf_normalized_sinc(2 * cutoff_frequency * (x - (filter_size - 1) / 2.0)) * lf_hamming_window(x, filter_size - 1);
    }
}

struct lf_filter
{
    LF_FLOATING_POINT *tail;

    LF_FLOATING_POINT *kernel_real;
    LF_FLOATING_POINT *kernel_imag;

    int impulse_response_length;
    int processing_chunk_length;

    LF_FLOATING_POINT *workspaces[4];
};

int lf_filter_conv_result_length(lf_filter *filter)
{
    return filter->impulse_response_length + filter->processing_chunk_length - 1;
}

lf_filter* lf_filter_init(const LF_FLOATING_POINT *impulse_response, int impulse_response_length, int processing_chunk_length)
{
    lf_filter *filter = (lf_filter*)malloc(sizeof(lf_filter));

    filter->impulse_response_length = impulse_response_length;
    filter->processing_chunk_length = processing_chunk_length;

    filter->tail = (LF_FLOATING_POINT*)calloc(impulse_response_length - 1, sizeof(LF_FLOATING_POINT));

    filter->kernel_real = (LF_FLOATING_POINT*)malloc(sizeof(LF_FLOATING_POINT) * lf_filter_conv_result_length(filter));
    filter->kernel_imag = (LF_FLOATING_POINT*)malloc(sizeof(LF_FLOATING_POINT) * lf_filter_conv_result_length(filter));
    lf_filter_replace_impulse_response(filter, impulse_response);
    
    for(int i = 0; i < 4; i++)
        filter->workspaces[i] = (LF_FLOATING_POINT*)malloc(sizeof(LF_FLOATING_POINT) * lf_filter_conv_result_length(filter));

    return filter;
}

void lf_filter_uninit(lf_filter* filter)
{
    free(filter->tail);
    free(filter->kernel_real);
    free(filter->kernel_imag);
    for(int i = 0; i < 4; i++)
        free(filter->workspaces[i]);
    free(filter);
}

void lf_filter_replace_impulse_response(lf_filter *filter, const LF_FLOATING_POINT *impulse_response)
{
    LF_FLOATING_POINT impulse_response_padded[lf_filter_conv_result_length(filter)];
    memset(impulse_response_padded, 0, sizeof(LF_FLOATING_POINT) * lf_filter_conv_result_length(filter));
    for(int i = 0; i < filter->impulse_response_length; i++)
        impulse_response_padded[i] = impulse_response[i], 0;

    LF_FLOATING_POINT all_zero_array[lf_filter_conv_result_length(filter)];
    memset(all_zero_array, 0, sizeof(LF_FLOATING_POINT) * lf_filter_conv_result_length(filter));

    lf_fft(impulse_response_padded, all_zero_array, filter->kernel_real, filter->kernel_imag, lf_filter_conv_result_length(filter));
}

void lf_filter_process(lf_filter* filter, const LF_FLOATING_POINT* input, LF_FLOATING_POINT* output)
{
    for(int i = 0; i < 4; i++)
        memset(filter->workspaces[i], 0, sizeof(LF_FLOATING_POINT) * lf_filter_conv_result_length(filter));

    for(int i = 0; i < filter->processing_chunk_length; i++)
        filter->workspaces[0][i] = input[i];

    lf_fft(filter->workspaces[0], filter->workspaces[1], filter->workspaces[2], filter->workspaces[3], lf_filter_conv_result_length(filter));

    for(int i = 0; i < lf_filter_conv_result_length(filter); i++)
        lf_complex_multiply(filter->kernel_real[i], filter->kernel_imag[i], filter->workspaces[2][i], filter->workspaces[3][i], filter->workspaces[0] + i, filter->workspaces[1] + i);

    lf_ifft(filter->workspaces[0], filter->workspaces[1], filter->workspaces[2], filter->workspaces[3], lf_filter_conv_result_length(filter));

    for(int i = 0; i < lf_filter_conv_result_length(filter); i++)
    {
        LF_FLOATING_POINT out = filter->workspaces[2][i];
        if(i < filter->impulse_response_length - 1)
            out += filter->tail[i];

        if(i < filter->processing_chunk_length)
            output[i] = out;
        if(i >= filter->processing_chunk_length)
            filter->tail[i - filter->processing_chunk_length] = out;
    }
}

void lf_windowed_sinc_lowpass(LF_FLOATING_POINT *output, int N, lf_window_func_t window_func, LF_FLOATING_POINT cutoff_freq)
{
    for(int i = 0; i < N; i++)
        output[i] = lf_windowed_sinc(i, cutoff_freq, window_func, N);
}

void lf_windowed_sinc_highpass(LF_FLOATING_POINT *output, int N, lf_window_func_t window_func, LF_FLOATING_POINT cutoff_freq)
{
    for(int i = 0; i < N; i++)
    {
        if(i % 2 == 0)
            output[i] = lf_windowed_sinc(i, 0.5 - cutoff_freq, window_func, N);
        else
            output[i] = -lf_windowed_sinc(i, 0.5 - cutoff_freq, window_func, N);
    }
}

void lf_windowed_sinc_bandpass(LF_FLOATING_POINT *output, int N, lf_window_func_t window_func, LF_FLOATING_POINT cutoff_freq, LF_FLOATING_POINT band_width)
{
    for(int i = 0; i < N; i++)
    {
        output[i] = lf_windowed_sinc(i, cutoff_freq + band_width / 2, window_func, N)
            - lf_windowed_sinc(i, cutoff_freq - band_width / 2, window_func, N);
    }
}

void lf_windowed_sinc_bandreject(LF_FLOATING_POINT *output, int N, lf_window_func_t window_func, LF_FLOATING_POINT cutoff_freq, LF_FLOATING_POINT band_width)
{
    for(int i = 0; i < N; i++)
    {
        if(i % 2 == 0)
            output[i] = lf_windowed_sinc(i, cutoff_freq - band_width / 2, window_func, N)
                + lf_windowed_sinc(i, 0.5 - (cutoff_freq + band_width / 2), window_func, N);
        else
            output[i] = lf_windowed_sinc(i, cutoff_freq - band_width / 2, window_func, N)
                - lf_windowed_sinc(i, 0.5 - (cutoff_freq + band_width / 2), window_func, N);
    }
}

void lf_moving_average(LF_FLOATING_POINT *output, int N)
{
    for(int i = 0; i < N; i++)
        output[i] = (LF_FLOATING_POINT)1 / N;
}
