#include <iostream>
#include <array>
#include <sstream>
#include <cstdio>
#include <cstring>

#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "./stb_image.h"
#include "./stb_image_write.h"

#define DJUTIL_NEEDS_linalg
#define DJUTIL_NEEDS_ezstr
#define DJUTIL_NEEDS_argparse

#include "./djutil.hpp"

namespace linalg = djutil::linalg;
namespace ezstr = djutil::ezstr;
namespace argparse = djutil::argparse;

using linalg::fvec3;

typedef long long lint;

inline lint imod(const lint& a, const lint& b){return ((a % b) + b) % b;}

lint iclamp(const lint& a, const lint& b){
    if (a < 0){return 0;}
    else if (a >= b){return b;}
    else {return a;}
}

inline float fclamp(const float& v, const float& l, const float& h) {
    return (v < h ? (v > l ? v : l) : h);
}

template <typename VT>
inline VT lerp(const VT& a, const VT& b, const VT& c) {
    return a+((b-a)*c);
}

struct heightmap_t {
    std::vector<float> points = {};
    lint width = 0, height = 0;

    fvec3 get(const lint& x, const lint& y) {
        //return this->points[(imod(y, this->height) * this->width) + imod(x, this->width)];
        //return this->points[(iclamp(y, this->height-1) * this->width) + iclamp(x, this->width-1)];
        const float
            pcx = lerp(-1.0f, 1.0f, (float(x)/(this->width-1))),
            pcy = lerp(-1.0f, 1.0f, (float(y)/(this->height-1))),
            pcz = this->points[(imod(y, this->height) * this->width) + imod(x, this->width)]
        ;
        return fvec3{pcx, pcy, pcz};
    }

    void load(const std::string& fn, const float heightmul=0.1f, const float zmin=-1.0f, const float zmax=1.0f){
        this->width = 0;
        this->height = 0;
        this->points.clear();

        const float hminv = 0.0f;
        const float hmaxv = float(USHRT_MAX);
        const float hdtv = hmaxv - hminv;

        int w, h, nc;
        uint16_t* px = stbi_load_16(fn.c_str(), &w, &h, &nc, 1);
        if (px == nullptr){
            std::stringstream errss;
            errss << "[FATAL ERROR]: Failed to load the input image at " << std::quoted(fn) << ": " << stbi_failure_reason() << "\n";
            throw std::runtime_error(errss.str());
        }
        this->width = w;
        this->height = h;
        this->points.resize(this->width * this->height);
        for (lint y = 0; y < this->height; y++){
            //const float pcy = lerp(-1.0f, 1.0f, (float(y)/(h-1)));
            for (lint x = 0; x < this->width; x++){
                //const float pcx = lerp(-1.0f, 1.0f, (float(x)/(w-1)));
                uint16_t inpx = px[(y * this->width) + x];
                const float pcz = lerp(zmin, zmax, ((float(inpx)-hminv)/hdtv)*heightmul);
                this->points[(y * this->width) + x] = pcz;
            }
        }

        STBI_FREE(px); px = nullptr;
    }

    heightmap_t() {}
    ~heightmap_t() {}
};

using rgb24_t = std::array<uint8_t, 3>;

void print_app_info_message() {
    std::cout << R"""(
    height2normal - a heightmap to normalmap converter.
    (c) 2023 rudolph4286, all lefts reserved.
    discord: rudolph4286
    email: fzerowipeoutlover1998@gmail.com OR djpaul3050@gmail.com
    github: github.com/SIGSEGV-666/

    This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
    You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.

)""";
}
void print_help_message() {
    std::cout << R"""(
    USAGE:
      height2normal INPUT_IMAGE OUTPUT_IMAGE [OPTIONS and/or FLAGS...]
    SYNTAX:
      height2normal INPUT_IMAGE OUTPUT_IMAGE -opt1 value1 --flags_can_be_here_too -opt2 value2 --flag2 -options_can_be_here_as_well value_is_here --flag3

    Converts a heightmap specified by the filename INPUT_IMAGE to a normalmap specified by the filename OUTPUT_IMAGE through the use of Sobel filtering.

    OPTIONS:
        -heightmul <int-or-decimal-value>:
            Heightmap value multiplier. Equal to 0.1 if omitted.
        -zmin <int-or-decimal-value>:
            Base minimum height for the heightmap. Equal to -1.0 if omitted.
        -zmax <int-or-decimal-value>:
            Base maximum height for the heightmap. Equal to 1.0 if omitted.
        -png_complevel <int-value>:
            Compression level to use if writing a .png file.
            Higher values mean more compression but takes longer to write the output PNG.
            Lower values mean less compression and shorter writing times, but the output PNG file will be larger.
            Equal to 8 if omitted.
    FLAGS:
        --xinv:
            Inverts the normal map's x axis (inverts the red channel values).
        --yinv:
            Inverts the normal map's y axis (inverts the green channel values).
        --zinv:
            Inverts the normal map's z axis (inverts the blue channel values).
        --tga_withrle:
            Uses RLE compression if writing a TGA file if specified.
        --help:
            Displays this text.
    Usage tips:
        By default, this tool generates OpenGL standard normalmaps.
        To generate a DirectX standard normalmap: pass just the --yinv flag to the program.
        You can alter the 'steepness' of the generated normal map by specifying a heightmul value greater than 0.1

        See the USAGE NOTES section of stb_image.h to see all of the image file formats this accepts as input.

        This tool can only write normal maps as either rgb24 PNGs or rgb24 TGAs
        (the 8-bits-per-channel limitation is stb_image_write.h's limitation, not mine.)

        The file extensions this program recognizes for output filenames are:
            '.png' or '.PNG' for PNG images
            '.tga' or '.TGA' for Truevision TGA images.

)""";
}

enum output_fmt : int {
   unknown = 0,
   png = 1,
   tga = 2
};
int main(int argc, const char* argv[]) {
    print_app_info_message();
    if (argc < 3){/*std::cerr << "Need at least 3 arguments!\n";*/ print_help_message(); return 0;}
    argparse::ArgumentParser<char> ap; ap.load_tokens(argv+3, argv+argc);
    heightmap_t hmap;
    std::string input_file = argv[1], output_file = argv[2];
    float heightmul = 0.1f;
    float xsign = 1.0f;
    float ysign = 1.0f;
    float zsign = 1.0f;

    float zmin = -1.0f;
    float zmax = 1.0f;

    int png_complevel = 8;
    bool tga_withrle = false;

    for (argparse::parsedarg_t<char> parg = {}; ap.adv(parg);){
        switch (parg.prefix){
            case argparse::argprefix::dash: {
                if (parg.name == "heightmul"){heightmul = std::stof(parg.value);}
                else if (parg.name == "zmin"){zmin = std::stof(parg.value);}
                else if (parg.name == "zmax"){zmax = std::stof(parg.value);}
                else if (parg.name == "png_complevel"){png_complevel = std::stoi(parg.value);}
                break;
            }
            case argparse::argprefix::double_dash: {
                if (parg.name == "xinv"){xsign = -1.0f;}
                else if (parg.name == "yinv"){ysign = -1.0f;}
                else if (parg.name == "zinv"){zsign = -1.0f;}
                else if (parg.name == "tga_withrle"){tga_withrle = true;}
                else if (parg.name == "help"){print_help_message(); return 0;}
                break;
            }
            default: ;
        }
    }

    stbi_write_png_compression_level = png_complevel;
    stbi_write_tga_with_rle = int(tga_withrle);

    std::cout << "[INFO]: Input (heightmap) filename = " << std::quoted(input_file) << "\n";
    std::cout << "[INFO]: Output (normalmap) filename = " << std::quoted(output_file) << "\n";

    output_fmt outfmt = output_fmt::unknown;
    std::string ofname, ofext; djutil::ezstr::rpartition<char>(ofname, ofext, output_file, ".");
    if (ofext == "png" || ofext == "PNG") {
        std::cout << "[INFO]: The normal map will be written as a PNG file.\n";
        outfmt = output_fmt::png;
    }
    else if (ofext == "tga" || ofext == "TGA"){
        std::cout << "[INFO]: The normal map will be written as a Truevision TGA file.\n";
        outfmt = output_fmt::tga;
    }
    else {
        std::cerr << "[FATAL ERROR]: Unknown normal map file extension: " << std::quoted(ofext) << ".\nPlease consult the --help dialog for a list of supported file extensions for normal maps.\n\n";
        return -1;
    }

    hmap.load(input_file, heightmul, zmin, zmax);
    std::cout << "[INFO]: Height map loaded successfully.\n";
    std::cout << "[INFO]: Image dimensions: w=" << hmap.width << ", h=" << hmap.height << "\n";
    std::vector<rgb24_t> normalmap = {}; normalmap.resize(hmap.width * hmap.height);

    const std::array<int, 2> sobel_offsets[9] = {
        {-1,1}, {0,1}, {1,1},
        {-1,0}, {0,0}, {1,0},
        {-1,-1}, {0,-1}, {1,-1}
    };
    float sobel_kernel[9] = {};

    std::cout << "[INFO]: Computing the normal map texels.\nPlease wait as this may take a few moments...\n";

    lint num_texels = hmap.width * hmap.height;
    size_t num_texels_strlen = std::to_string(num_texels).size();
    lint curtexel = 0;
    lint progress = 0, last_progress = 0;
    std::cout << "have to compute " << num_texels << " texels.\n\n";
    //std::cout << "Progress:   0 %";
    size_t progress_str_length = 0;
    for (lint y = 0; y < hmap.height; y++){
        for (lint x = 0; x < hmap.width; x++){
            last_progress = progress;
            curtexel++;
            progress = (curtexel * 100) / num_texels;
            if (progress > last_progress){
                for (int k = 0; k < progress_str_length; k++){std::cout << "\b \b";}
                //memset(percentage_str+2, 0, 6);
                std::stringstream pss;
                int pc = int(progress/5);

                pss << "Progress: " << std::setw(num_texels_strlen) << curtexel << "/" << num_texels << " texels [" ;
                for (int p = 0; p < 20; p++){
                    pss << (p <= pc ? '=' : ' ');
                }
                pss << "] " << std::setw(3) << progress << "%";
                std::cout << pss.str();
                std::cout.flush();
                progress_str_length = pss.str().size();

            }
            for (int i = 0; i < 9; i++){
                const auto& ofs = sobel_offsets[i];
                sobel_kernel[i] = hmap.get(x+ofs[0], y+ofs[1])[2];
            }

            const fvec3 normal = linalg::normalize(fvec3{
                60.0f * -(sobel_kernel[2]-sobel_kernel[0]+2*(sobel_kernel[5]-sobel_kernel[3])+sobel_kernel[8]-sobel_kernel[6]),
                60.0f * -(sobel_kernel[6]-sobel_kernel[0]+2*(sobel_kernel[7]-sobel_kernel[1])+sobel_kernel[8]-sobel_kernel[2]),
                1.0f
            });
            rgb24_t& outn = normalmap[(y * hmap.width) + x];
            outn[0] = uint8_t(fclamp(lerp(127.0f, 255.0f, normal[0] * xsign), 0.0f, 255.0f));
            outn[1] = uint8_t(fclamp(lerp(127.0f, 255.0f, normal[1] * ysign), 0.0f, 255.0f));
            outn[2] = uint8_t(fclamp(lerp(127.0f, 255.0f, normal[2] * zsign), 0.0f, 255.0f));
        }
    }
    std::cout << "\n\n";
    std::cout << "[INFO]: Normal map generation finished!\n";
    switch (outfmt){
        case output_fmt::png:
            std::cout << "[INFO]: Writing the normal map as a PNG file.\n";
            stbi_write_png(output_file.c_str(), hmap.width, hmap.height, 3, normalmap.data(), hmap.width * 3);
            break;
        case output_fmt::tga:
            std::cout << "[INFO]: Writing the normal map as a TGA file.\n";
            stbi_write_tga(output_file.c_str(), hmap.width, hmap.height, 3, normalmap.data());
            break;
        default: ;
    }
    std::cout << "Done!\n";
    return 0;
}

