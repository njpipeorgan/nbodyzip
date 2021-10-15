// nbodyzip, a lossy compressor of dark matter-only N-body simulation data
//
// Copyright 2021 Tianhuan Lu
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//

#include <cmath>
#include <cstdio>
#include <array>
#include <random>
#include <tuple>
#include <vector>

void _swap_bytes_2(uint8_t* data) {
    std::swap(data[0], data[1]);
}

void _swap_bytes_4(uint8_t* data) {
    std::swap(data[0], data[3]);
    std::swap(data[1], data[2]);
}

void _swap_bytes_8(uint8_t* data) {
    std::swap(data[0], data[7]);
    std::swap(data[1], data[6]);
    std::swap(data[2], data[5]);
    std::swap(data[3], data[4]);
}

template<typename T>
void swap_bytes(T* data) {
    if constexpr (sizeof(T) == 2) {
        _swap_bytes_2(reinterpret_cast<uint8_t*>(data));
    } else if constexpr (sizeof(T) == 4) {
        _swap_bytes_4(reinterpret_cast<uint8_t*>(data));
    } else if constexpr (sizeof(T) == 8) {
        _swap_bytes_8(reinterpret_cast<uint8_t*>(data));
    } else {
        static_assert(sizeof(T) == 2, "cannot do byte swap");
    }
}

struct bit_stream {
    uint8_t* init_{};
    uint8_t* begin_{};
    uint8_t* end_{};
    size_t cur_bit_{};

    bit_stream(void* begin, void* end) :
        init_{(uint8_t*)begin}, begin_{(uint8_t*)begin}, end_{(uint8_t*)end} {
    }

    void reset(void* begin, void* end) {
        begin_ = (uint8_t*)begin + (begin_ - init_);
        end_ = (uint8_t*)end;
        init_ = (uint8_t*)begin;
    }

    size_t tell() const {
        return 8u * size_t(begin_ - init_) + cur_bit_;
    }

    uint64_t read(size_t n) {
        auto new_bit = cur_bit_ + n % 8u;
        auto new_begin = begin_ + n / 8u + new_bit / 8u;
        new_bit %= 8u;
        if (n == 0u || n > 64u) {
            fprintf(stderr, "error: cannot read 0 bits or more than 64 bits\n");
            std::exit(1);
        }
        if (new_begin >= end_ && (new_begin > end_ || (new_begin == end_ && new_bit > 0u))) {
            fprintf(stderr, "error: bit stream out of range.\n");
            std::exit(1);
        }
        auto ret = uint64_t{};
        for (; begin_ < new_begin; ++begin_) {
            ret = (ret >> 8u) | (uint64_t(*begin_) << 56u);
        }
        if (new_bit > 0u) {
            ret = (ret >> new_bit) | (uint64_t(*begin_) << (64u - new_bit));
        }
        ret = ret >> (64u - n);
        cur_bit_ = new_bit;
        return ret;
    }

    uint64_t read() {
        if (begin_ >= end_) {
            fprintf(stderr, "error: bit stream out of range.\n");
            std::exit(1);
        }
        if (cur_bit_ < 7) {
            auto ret = (uint64_t(*begin_) >> cur_bit_) & uint64_t(1);
            ++cur_bit_;
            return ret;
        } else {
            auto ret = (uint64_t(*begin_) >> 7u) & uint64_t(1);
            ++begin_;
            cur_bit_ = 0;
            return ret;
        }
    }

    void write(size_t n, uint64_t data) {
        auto new_bit = cur_bit_ + n % 8u;
        auto new_begin = begin_ + n / 8u + new_bit / 8u;
        new_bit %= 8u;
        if (n == 0u || n > 64u) {
            fprintf(stderr, "error: cannot write 0 bits or more than 64 bits\n");
            std::exit(1);
        }
        if (new_begin >= end_ && (new_begin > end_ || (new_begin == end_ && new_bit > 0u))) {
            fprintf(stderr, "error: bit stream out of range.\n");
            std::exit(1);
        }
        for (; begin_ < new_begin; ++begin_, cur_bit_ = 0) {
            *begin_ |= (data << (56u + cur_bit_)) >> 56u;
            data >>= (8u - cur_bit_);
        }
        if (new_bit > 0u) {
            *begin_ |= (data << (64u - new_bit + cur_bit_)) >> (64u - new_bit);
        }
        cur_bit_ = new_bit;
    }

    void write(uint64_t data) {
        *begin_ |= (data & uint64_t(1)) << cur_bit_;
        if (cur_bit_ < 7) {
            ++cur_bit_;
        } else {
            cur_bit_ = 0;
            ++begin_;
        }
    }
};

struct posvel_data_t {
    double time{};
    uint64_t num_particles{};
    double mass{};
    double eps{};
    std::vector<std::array<uint32_t, 3>> pos_data{};
    std::vector<std::array<float, 3>> vel_data{};

    void read(FILE* str) {
        fprintf(stderr, "info: reading POSVEL file...\n");
        std::fread(this, 8u, 4u, str);
        pos_data.resize(num_particles);
        vel_data.resize(num_particles);
        std::fread(pos_data.data(), sizeof(pos_data[0]), num_particles, str);
        std::fread(vel_data.data(), sizeof(vel_data[0]), num_particles, str);
        fprintf(stderr, "info: time = %g, N = %d, mass = %g, eps = %g\n", time, int(num_particles), mass, eps);
    }

    void write(FILE* str) const {
        fprintf(stderr, "info: writing POSVEL file...\n");
        std::fwrite(this, 8u, 4u, str);
        std::fwrite(pos_data.data(), sizeof(pos_data[0]), num_particles, str);
        std::fwrite(vel_data.data(), sizeof(vel_data[0]), num_particles, str);
    }
};

struct zip_data_t {
    double time;
    uint64_t num_particles{};
    double mass{};
    double eps{};
    double vel_scale{};
    uint64_t num_bins{};
    uint64_t zip_counts_size{}; // size in bytes
    std::vector<uint32_t> counts{};
    std::vector<uint64_t> zip_counts{};
    std::vector<uint32_t> pos_data; // x: 11 bits, y: 11 bits, z: 10 bits
    std::vector<std::array<int16_t, 3>> vel_data;

    static constexpr double ARCTAN_FACTOR = 1.571;
    static constexpr double VEL_ENCODE_FACTOR = 32768.0;

    void read(FILE* str) {
        fprintf(stderr, "info: reading ZIP file...\n");
        std::fread(&time, 8u, 1u, str);
        std::fread(&num_particles, 8u, 1u, str);
        std::fread(&mass, 8u, 1u, str);
        std::fread(&eps, 8u, 1u, str);
        std::fread(&vel_scale, 8u, 1u, str);
        std::fread(&num_bins, 8u, 1u, str);
        std::fread(&zip_counts_size, 8u, 1u, str);

        if (num_bins * num_bins * num_bins != num_particles) {
            fprintf(stderr, "error: the number of particles does not match the number of bins.\n");
            std::exit(1);
        }

        zip_counts.resize(zip_counts_size / sizeof(zip_counts[0]) + 1);
        pos_data.resize(num_particles);
        vel_data.resize(num_particles);
        std::fread(zip_counts.data(), 1u, zip_counts_size, str);
        std::fread(pos_data.data(), sizeof(pos_data[0]), num_particles, str);
        std::fread(vel_data.data(), sizeof(vel_data[0]), num_particles, str);

        // decompress counts
        // 0         : 0
        // 1         : 10
        // 2         : 110
        // 3         : 1110
        // 4-7       : 1111 0xx
        // 8-15      : 1111 10xx x
        // 16-31     : 1111 110x xxx
        // 32-255    : 1111 1110 xxxx xxxx
        // 256-65535 : 1111 1111 0xxx xxxx xxxx xxxx x
        // 65536-    : 1111 1111 10xx xxxx xxxx xxxx xxxx xxxx xxxx xxxx xx

        counts.clear();
        counts.resize(num_particles);
        auto bit_str = bit_stream(zip_counts.data(), zip_counts.data() + zip_counts.size());
        for (size_t i = 0; i < num_particles; ++i) {
            size_t ones = 0;
            for (; bit_str.read(); ++ones) {}
            if (ones < 4u) {
                counts[i] = uint32_t(ones);
            } else if (ones < 10u) {
                counts[i] = uint32_t(
                    ones == 4u ? 4u +  bit_str.read(2)  :
                    ones == 5u ? 8u +  bit_str.read(3)  :
                    ones == 6u ? 16u + bit_str.read(4)  :
                    ones == 7u ? 0u +  bit_str.read(8)  :
                    ones == 8u ? 0u +  bit_str.read(16) :
                                 0u +  bit_str.read(32)
                );
            } else {
                fprintf(stderr, "error: invalid compress count encoding.\n");
                std::exit(1);
            }
        }
    }

    void write(FILE* str) {
        fprintf(stderr, "info: writing ZIP file...\n");
        std::fwrite(&time, 8u, 1u, str);
        std::fwrite(&num_particles, 8u, 1u, str);
        std::fwrite(&mass, 8u, 1u, str);
        std::fwrite(&eps, 8u, 1u, str);
        std::fwrite(&vel_scale, 8u, 1u, str);
        std::fwrite(&num_bins, 8u, 1u, str);

        // compress counts
        zip_counts.clear();
        zip_counts.resize(std::max(size_t(1000), num_particles / 2u));
        auto bit_str = bit_stream(zip_counts.data(), zip_counts.data() + zip_counts.size());
        for (size_t i = 0; i < num_particles; ++i) {
            const auto count = counts[i];
            if (count == 0u) {
                bit_str.write(0b0u);
            } else if (count == 1u) {
                bit_str.write(2u, 0b01ull);
            } else if (count == 2u) {
                bit_str.write(3u, 0b011ull);
            } else if (count == 3u) {
                bit_str.write(4u, 0b0111ull);
            } else if (count < 8u) {
                bit_str.write(7u, 0b0111'1ull + ((count - 4u) << 5u));
            } else if (count < 16u) {
                bit_str.write(9u, 0b0111'11ull + ((count - 8u) << 6u));
            } else if (count < 32u) {
                bit_str.write(11u, 0b0111'111ull + ((count - 16u) << 7u));
            } else if (count < 256u) {
                bit_str.write(16u, 0b0111'1111ull + (count << 8u));
            } else if (count < 65536u) {
                bit_str.write(25u, 0b0111'1111'1ull + (count << 9u));
            } else {
                bit_str.write(42u, 0b0111'1111'11ull + (count << 10u));
            }
        }
        const auto zip_counts_size = bit_str.tell() / 8u + 1u;
        zip_counts.resize((zip_counts_size - 1u) / 8u + 1u);
        fprintf(stderr, "info: compress count, n_bytes = %g, bits/particle = %g\n",
            double(zip_counts_size), double(zip_counts_size * 8.0 / num_particles));

        std::fwrite(&zip_counts_size, 8u, 1u, str);
        std::fwrite(zip_counts.data(), 1u, zip_counts_size, str);
        std::fwrite(pos_data.data(), sizeof(pos_data[0]), num_particles, str);
        std::fwrite(vel_data.data(), sizeof(vel_data[0]), num_particles, str);
    }

    void level_check(uint64_t level) const {
        if (num_bins * num_bins * num_bins != num_particles) {
            fprintf(stderr, "error: the number of particles does not match the number of bins.\n");
            fprintf(stderr, "       num_particles = %ld, num_bins = %ld\n", num_particles, num_bins);
            std::exit(1);
        }
        if ((1ull << level) != num_bins || level > 10u) {
            fprintf(stderr, "error: the number of bins is not a power (<= 10) of 2.\n");
            fprintf(stderr, "       num_particles = %ld, num_bins = %ld\n", num_particles, num_bins);
            std::exit(1);
        }
    }

    void to_posvel(posvel_data_t& posvel) {
        posvel.time = time;
        posvel.num_particles = num_particles;
        posvel.mass = mass;
        posvel.eps = eps;
        const auto level = (uint64_t)std::llround(std::log2(double(num_bins)));
        level_check(level);

        posvel.pos_data.resize(num_particles);
        posvel.vel_data.resize(num_particles);
        auto i = size_t{};
        auto partial_sum = 0;
        for (size_t bin = 0; bin < num_particles; ++bin) {
            partial_sum += counts[bin];
        }

        for (size_t bin = 0; bin < num_particles; ++bin) {
            auto i_end = i + size_t(counts[bin]);
            for (; i < i_end; ++i) {
                const auto pos_code = pos_data[i];
                const uint32_t bin_x = uint32_t(bin) & ((1u << level) - 1u);
                const uint32_t bin_y = (uint32_t(bin) >> level) & ((1u << level) - 1u);
                const uint32_t bin_z = uint32_t(bin) >> (level + level);
                const uint32_t code_x = pos_code & 0b0111'1111'1111u;
                const uint32_t code_y = (pos_code >> 11u) & 0b0111'1111'1111u;
                const uint32_t code_z = pos_code >> 22u;
                const uint32_t delta = 0x7fff'ffffu >> (level + 10u);

                auto& p_pos = posvel.pos_data[i];
                auto& p_vel = posvel.vel_data[i];
                p_pos[0] = (bin_x << (32u - level)) | (code_x << (32u - 11u - level)) | (delta >> 1u);
                p_pos[1] = (bin_y << (32u - level)) | (code_y << (32u - 11u - level)) | (delta >> 1u);
                p_pos[2] = (bin_z << (32u - level)) | (code_z << (32u - 10u - level)) | (delta);
                for (int iv = 0; iv < 3; ++iv) {
                    p_vel[iv] = (vel_scale / ARCTAN_FACTOR) *
                        std::tan(vel_data[i][iv] / (VEL_ENCODE_FACTOR / ARCTAN_FACTOR));
                }
            }
        }
    }

    void from_posvel(const posvel_data_t& posvel) {
        time = posvel.time;
        num_particles = posvel.num_particles;
        mass = posvel.mass;
        eps = posvel.eps;
        num_bins = (uint64_t)std::llround(std::pow(double(num_particles), 1.0 / 3.0));
        const auto level = (uint64_t)std::llround(std::log2(double(num_bins)));
        level_check(level);

        counts.clear();
        counts.resize(num_particles);
        pos_data.resize(num_particles);
        vel_data.resize(num_particles);
        auto which_bin = std::vector<uint32_t>(num_particles);
        auto sum_counts = std::vector<uint32_t>(num_particles);
        auto sum_v2 = 0.0;

        for (size_t i = 0; i < num_particles; ++i) {
            auto& p_pos = posvel.pos_data[i];
            auto& p_vel = posvel.vel_data[i];
            size_t bin = (p_pos[0] >> (32u - level)) |
                ((p_pos[1] >> (32u - level)) << level) |
                ((p_pos[2] >> (32u - level)) << (level + level));
            which_bin[i] = uint32_t(bin);
            counts[bin] += 1u;
            sum_v2 += p_vel[0] * p_vel[0] + p_vel[1] * p_vel[1] + p_vel[2] * p_vel[2];
        }
        for (size_t i = 1; i < num_particles; ++i) {
            sum_counts[i] = sum_counts[i - 1] + counts[i - 1];
        }
        vel_scale = std::sqrt(sum_v2 / (3.0 * num_particles));
        const auto inv_vel_scale = 1.0 / vel_scale;

        for (size_t i = 0; i < num_particles; ++i) {
            auto& p_pos = posvel.pos_data[i];
            auto& p_vel = posvel.vel_data[i];
            auto pos_code =
                (((p_pos[0] >> (32u - 11u - level)) & uint32_t(0b0111'1111'1111)) << 0u ) |
                (((p_pos[1] >> (32u - 11u - level)) & uint32_t(0b0111'1111'1111)) << 11u) |
                (((p_pos[2] >> (32u - 10u - level)) & uint32_t(0b0011'1111'1111)) << 22u);
            auto& insert_pos = sum_counts[which_bin[i]];
            pos_data[insert_pos] = pos_code;
            for (size_t iv = 0; iv < 3; ++iv) {
                vel_data[insert_pos][iv] = static_cast<int16_t>(std::lround(
                    (VEL_ENCODE_FACTOR / ARCTAN_FACTOR) *
                    std::atan(p_vel[iv] * inv_vel_scale * ARCTAN_FACTOR)
                ));
            }
            ++insert_pos;
        }
    }
};

struct tipsy_data_t {
    struct particle_t {
        float mass{};
        std::array<float, 3> pos{};
        std::array<float, 3> vel{};
        float eps{};
        float phi{};
    };

    double time;
    uint32_t num_particles;
    std::vector<particle_t> particles{};

    static void swap_particle_bytes(particle_t& p) {
        swap_bytes(&p.mass);
        swap_bytes(&p.pos[0]);
        swap_bytes(&p.pos[1]);
        swap_bytes(&p.pos[2]);
        swap_bytes(&p.vel[0]);
        swap_bytes(&p.vel[1]);
        swap_bytes(&p.vel[2]);
        swap_bytes(&p.eps);
        swap_bytes(&p.phi);
    }

    void read(FILE* str) {
        fprintf(stderr, "info: reading TIPSY file...\n");
        std::fread(&time, 8u, 1u, str);
        std::fread(&num_particles, 4u, 1u, str);
        std::fseek(str, 20, SEEK_CUR);
        swap_bytes(&time);
        swap_bytes(&num_particles);
        particles.resize(num_particles);
        std::fread(particles.data(), sizeof(particle_t), num_particles, str);
        for (auto& p : particles) {
            swap_particle_bytes(p);
        }
        fprintf(stderr, "info: time = %g, N = %d, mass = %g, eps = %g\n",
            time, int(num_particles), particles.at(0).mass, particles.at(0).eps);
    }

    void write(FILE* str) const {
        fprintf(stderr, "info: writing TIPSY file...\n");
        auto xdr_time = this->time;
        auto xdr_num_particles = this->num_particles;
        auto xdr_num_dimensions = uint32_t{3};
        const auto padding = uint32_t{};
        swap_bytes(&xdr_time);
        swap_bytes(&xdr_num_particles);
        swap_bytes(&xdr_num_dimensions);
        std::fwrite(&xdr_time, 8u, 1u, str);
        std::fwrite(&xdr_num_particles, 4u, 1u, str);
        std::fwrite(&xdr_num_dimensions, 4u, 1u, str);
        std::fwrite(&padding, 4u, 1u, str);
        std::fwrite(&xdr_num_particles, 4u, 1u, str);
        std::fwrite(&padding, 4u, 1u, str);
        std::fwrite(&padding, 4u, 1u, str);

        constexpr auto block_size = size_t(100);
        auto xdr_particles = std::vector<particle_t>(block_size);
        auto begin = particles.begin();
        while (begin < particles.end()) {
            auto end = std::min(begin + block_size, particles.end());
            std::copy(begin, end, xdr_particles.begin());
            for (auto& p : xdr_particles) {
                swap_particle_bytes(p);
            }
            std::fwrite(xdr_particles.data(), sizeof(particle_t), end - begin, str);
            begin = end;
        }
    }

    auto to_std_pos(const std::array<float, 3>& pos) {
        auto ret = std::array<uint32_t, 3>{};
        for (size_t i = 0; i < 3; ++i) {
            const auto intpos = std::llroundf(pos[i] * 4294967296.0f);
            ret[i] = static_cast<uint32_t>(intpos + 2147483648ll);
        }
        return ret;
    }

    auto from_std_pos(const std::array<uint32_t, 3>& pos) {
        auto ret = std::array<float, 3>{};
        for (size_t i = 0; i < 3; ++i) {
            const auto intpos = (static_cast<long long>(pos[i]) - 2147483648ll);
            ret[i] = static_cast<float>(intpos) / 4294967296.0f;
        }
        return ret;
    }

    void to_posvel(posvel_data_t& posvel) {
        posvel.time = this->time;
        posvel.num_particles = this->num_particles;
        posvel.mass = this->particles.at(0).mass;
        posvel.eps = this->particles.at(0).eps;
        posvel.pos_data.resize(num_particles);
        posvel.vel_data.resize(num_particles);
        for (size_t i = 0; i < num_particles; ++i) {
            posvel.pos_data[i] = to_std_pos(this->particles[i].pos);
            posvel.vel_data[i] = this->particles[i].vel;
        }
    }

    void from_posvel(const posvel_data_t& posvel) {
        this->time = posvel.time;
        this->num_particles = posvel.num_particles;
        this->particles.resize(num_particles);
        for (size_t i = 0; i < num_particles; ++i) {
            this->particles[i].mass = posvel.mass;
            this->particles[i].pos = from_std_pos(posvel.pos_data[i]);
            this->particles[i].vel = posvel.vel_data[i];
            this->particles[i].eps = posvel.eps;
            this->particles[i].phi = 0.0f;
        }
    }
};

auto random_id() {
    static std::random_device rd;
    static auto eng = std::mt19937_64{rd()};
    auto dist = std::uniform_int_distribution<>(int('A'), int('Z'));
    
    auto id = std::string(16, ' ');
    for (auto& ch : id) {
        ch = char(dist(eng));
    }
    return id;
}

enum class format_t {
    TIPSY, POSVEL, ZIP, INVALID
};

auto str_to_format(const std::string& arg) {
    return (
        arg == "tipsy"  || arg == "TIPSY" ? format_t::TIPSY  :
        arg == "posvel" || arg == "POVEL" ? format_t::POSVEL :
        arg == "zip"    || arg == "ZIP"   ? format_t::ZIP    :
        format_t::INVALID
    );
}

int main(int argc, char* argv[]) {

    const auto input_format = str_to_format(argc == 5 ? argv[1] : "");
    const auto output_format = str_to_format(argc == 5 ? argv[2] : "");
    const auto input_filename = std::string(argc == 5 ? argv[3] : "");
    auto output_filename = std::string(argc == 5 ? argv[4] : "");

    if (input_format == format_t::INVALID || output_format == format_t::INVALID) {
        printf(R"(nbodyzip, a lossy compressor of dark matter-only N-body simulation data
Usage:
    nbodyzip [input_format] [output_format] [input_file] [output_file]
input/output format can be:
    tipsy   (dark matter only)
    posvel  binary file storing positions and velocities
    zip     file compressed by nbodyzip
)");
        return 0;
    }

    auto input_stream = std::fopen(input_filename.c_str(), "rb");
    if (!input_stream) {
        fprintf(stderr, "error: failed to open file, %s\n", input_filename.c_str());
        return 1;
    }
    {
        auto output_stream = std::fopen(output_filename.c_str(), "rb");
        if (output_stream) {
            fprintf(stderr, "warning: output file already exists, %s\n", output_filename.c_str());
            output_filename += "." + random_id();
            fprintf(stderr, "         writing to a new file, %s\n", output_filename.c_str());
        }
    }
    auto output_stream = std::fopen(output_filename.c_str(), "wb");
    if (!output_stream) {
        fprintf(stderr, "error: failed to open file %s\n", output_filename.c_str());
        return 1;
    }

    auto posvel_buffer = posvel_data_t{};
    auto tipsy_buffer = tipsy_data_t{};
    auto zip_buffer = zip_data_t{};

    if (input_format == format_t::TIPSY) {
        tipsy_buffer.read(input_stream);
    } else if (input_format == format_t::POSVEL) {
        posvel_buffer.read(input_stream);
    } else if (input_format == format_t::ZIP) {
        zip_buffer.read(input_stream);
    }
    if (std::feof(input_stream)) {
        fprintf(stderr, "error: failed to read the input file.\n");
        return 1;
    }
    if (std::fseek(input_stream, 1, SEEK_CUR)) {
        fprintf(stderr, "warning: the end of the file is not reached at the end of reading.\n");
    }
    std::fclose(input_stream);

    if (input_format == format_t::TIPSY) {
        tipsy_buffer.to_posvel(posvel_buffer);
    } else if (input_format == format_t::ZIP) {
        zip_buffer.to_posvel(posvel_buffer);
    }

    if (output_format == format_t::TIPSY) {
        tipsy_buffer.from_posvel(posvel_buffer);
        tipsy_buffer.write(output_stream);
    } else if (output_format == format_t::POSVEL) {
        posvel_buffer.write(output_stream);
    } else if (output_format == format_t::ZIP) {
        zip_buffer.from_posvel(posvel_buffer);
        zip_buffer.write(output_stream);
    }
    if (std::ferror(output_stream)) {
        fprintf(stderr, "error: failed to write the file.\n");
        return 1;
    }
    std::fclose(output_stream);

    return 0;
}

