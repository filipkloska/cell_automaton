#include <curand_kernel.h>
typedef unsigned char u8;
#define EMPTY 0
#define OCEAN 1
#define SEA 2
#define COAST 3
#define LOWLANDS 4
#define HIGHLANDS 5
#define MOUNTAINS 6
#define CHECKED 7


// Struktura odpowiadająca GpuCoord z Rusta
struct GpuCoord {
    int y;
    int x;
};

extern "C" __global__ void setup_rand(
    curandState *state, 
    unsigned int seed, 
    int width, 
    int height) 
{
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    if (id < width * height) {
        curand_init(seed, id, 0, &state[id]);
    }
}


__constant__ double GPU_CHANCES[6][6] = {
    {0.8, 0.2, 0.1, 0.0, -0.1, -0.1}, 
    {0.2, 0.2, 0.8, 0.0, -0.1, -0.1},
    {0.1, 0.2, 0.1, 0.6, 0.0, 0.0},
    {0.0, 0.0, 0.1, 0.5, 0.2, 0.1},
    {-0.1, -0.1, 0.0, 0.4, 0.4, 0.2},
    {-0.1, -0.1, 0.0, 0.1, 0.7, 0.2}
};


extern "C" __global__ void biome_step(
    const u8* input,
    u8* output,
    curandState* states,
    int width,
    int height,
    int read_window_height,
    int read_window_width,
    int write_y_start, 
    int write_y_end) 
{
    //coords of the pixel
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;

    //sanity check
    if (x >= width || y >= height)
    {
        return;
    }
    //get the value from the coords + check if it's already filled
    int idx = y * width + x;
    u8 current = input[idx];
    if (current > 0 && current <= 6) {
        output[idx] = current;
        return;
    }
    //avoid safety offset areas
    if (y < write_y_start || y >= write_y_end ) {
        output[idx] = current; 
        return;
    }
    
    //sanity check 2 - we need at least one immediate neighbor to grow
    bool has_immediate = false;
    if (y > 0 && input[(y-1) * width + x] > 0)
    {
        has_immediate = true;
    }         
    else if (y < height - 1 && input[(y+1) * width + x] > 0)
    {
        has_immediate = true;
    } 
    else if (x > 0 && input[y * width + (x-1)] > 0)
    {
        has_immediate = true;
    }
    else if (x < width - 1 && input[y * width + (x+1)] > 0)
    {
        has_immediate = true;
    }
    

    if (!has_immediate) {
        output[idx] = 0;
        return;
    }

    //get the accumulated chances from neighbors
    float chances[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    for (int wy = -read_window_height; wy <= read_window_height; wy++) {
        for (int wx = -read_window_width; wx <= read_window_width; wx++) {
            if (wx == 0 && wy == 0)
            {
                continue;
            } 
            int nx = x + wx;
            int ny = y + wy;

            if (nx >= 0 && nx < width && ny >= 0 && ny < height) {
                u8 neighbor_val = input[ny * width + nx];
                if (neighbor_val >= 1 && neighbor_val <= 6) {
                    int neighbor_idx = neighbor_val - 1;
                    //further neighbors have less impact
                    double mult = (abs(wx) + abs(wy) > 1) ? 0.25 : 1.0; 

                    for (int b = 0; b < 6; b++) {
                        chances[b] += GPU_CHANCES[neighbor_idx][b] * mult;
                    }
                }
            }
        }
    }

    float sum = 0;
    for (int i = 0; i < 6; i++) {
        if (chances[i] < 0.0) 
        {
            chances[i] = 0.0;
        }
        sum += chances[i];
    }

    curandState localState = states[idx];
    float rand_val = curand_uniform(&localState);
    states[idx] = localState;
    u8 selected = 0;
    
    if (sum <= 0.0) {
        float running_sum = 0;
        for (int i = 0; i < 6; i++) {
            running_sum += (1.0 / 6.0);
            if (rand_val <= running_sum) {
                selected = (u8)(i + 1);
                break;
            }
        }
    }
    else
    {
        float running_sum = 0;
        for (int i = 0; i < 6; i++) {
            running_sum += (chances[i] / sum);
            if (rand_val <= running_sum) {
                selected = (u8)(i + 1);
                break;
            }
        }
    }
    output[idx] = selected;
}
