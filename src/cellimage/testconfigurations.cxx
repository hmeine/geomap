#include <iostream>
#include <string>
#include "celltypes.hxx"

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
