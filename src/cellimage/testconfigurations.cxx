/************************************************************************/
/*                                                                      */
/*               Copyright 2007-2019 by Hans Meine                      */
/*                                                                      */
/*    Permission is hereby granted, free of charge, to any person       */
/*    obtaining a copy of this software and associated documentation    */
/*    files (the "Software"), to deal in the Software without           */
/*    restriction, including without limitation the rights to use,      */
/*    copy, modify, merge, publish, distribute, sublicense, and/or      */
/*    sell copies of the Software, and to permit persons to whom the    */
/*    Software is furnished to do so, subject to the following          */
/*    conditions:                                                       */
/*                                                                      */
/*    The above copyright notice and this permission notice shall be    */
/*    included in all copies or substantial portions of the             */
/*    Software.                                                         */
/*                                                                      */
/*    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND    */
/*    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   */
/*    OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          */
/*    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT       */
/*    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,      */
/*    WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING      */
/*    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR     */
/*    OTHER DEALINGS IN THE SOFTWARE.                                   */
/*                                                                      */
/************************************************************************/

#include <iostream>
#include <string>
#include "cellconfigurations.hxx"

int main()
{
    std::cerr << "Errors: " << std::endl;

    for(int configuration = 0; configuration<256; ++configuration)
    {
        int bits[8];
        for(int bitpos=0; bitpos<8; ++bitpos)
            bits[bitpos] = (configuration >> bitpos) & 1;

        int reflected  =  bits[0] |
						 (bits[7] << 1) |
						 (bits[6] << 2) |
						 (bits[5] << 3) |
						 (bits[4] << 4) |
						 (bits[3] << 5) |
						 (bits[2] << 6) |
						 (bits[1] << 7);

        int other = configuration;

        for(int k=0; k<4; ++k)
        {
            if(cellConfigurations[configuration] != cellConfigurations[other])
            {
                std::cerr << configuration << ' ' << other << std::endl;
            }
            if(cellConfigurations[configuration] != cellConfigurations[reflected])
            {
                std::cerr << configuration << ' ' << reflected << std::endl;
            }

            other = ((other << 2) | (other >> 6)) & 0xff;

            reflected = ((reflected << 2) | (reflected >> 6)) & 0xff;
        }
    }
}
