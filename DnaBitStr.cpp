//	Metal - A fast methylation alignment and calling tool for WGBS data.
//	Copyright (C) 2017  Jonas Fischer
//
//	This program is free software: you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation, either version 3 of the License, or
//	(at your option) any later version.
//
//	This program is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//	GNU General Public License for more details.
//
//	You should have received a copy of the GNU General Public License
//	along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//	Jonas Fischer	jonaspost@web.de

#include "DnaBitStr.h"


DnaBitStr::DnaBitStr(const unsigned int s) :
        size(s)
    ,   bitSeq((size / 32) + 1)
    ,   bitMask((size / 32) + 1)
    ,   bitRevMask((size / 32) + 1)
{
}



void DnaBitStr::setBitStrLast(std::string&& seq)
{


    uint64_t bitStr = 0;
    uint64_t bitM = 0xffffffffffffffffULL;
    uint64_t bitRevM = 0xffffffffffffffffULL;
    for (size_t i = 1; i <= seq.size(); ++i)
    {

        const unsigned int shift = 64 - 2*i;
        // we don't need A here
        switch (seq[i - 1]){

            case 'C':
                {
                    bitStr |= (1ULL << shift);
                    bitM ^= (2ULL << shift);
                }
                break;

            case 'G':
                {
                    bitStr |= (2ULL << shift);
                    bitRevM ^= (2ULL << shift);
                }
                break;

            case 'T':
                bitStr |= (3ULL << (64 - 2*i));
                break;

            // unknowns, A and N will be encoded as zero; nothing to do
            default:
                break;
        }
    }
    bitSeq[size/32] = bitStr;
    bitMask[size/32] = bitM;
    bitRevMask[size/32] = bitRevM;
}
