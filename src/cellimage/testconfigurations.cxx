#include <iostream.h>
#include <string>
#include "configurations.h"

int main()
{
    int configuration = 0;

    cerr << "Test: " << cellConfigurations[0] << ' ' <<
                        cellConfigurations[255] << endl;

    for(; configuration<256; ++configuration)
    {
        int k;

        int bits[8];
        for(k=0; k<8; ++k)
        {
            bits[k] = (configuration >> k) & 1;
        }

        int reflected  =  bits[0] |
						 (bits[7] << 1) |
						 (bits[6] << 2) |
						 (bits[5] << 3) |
						 (bits[4] << 4) |
						 (bits[3] << 5) |
						 (bits[2] << 6) |
						 (bits[1] << 7);

        int r = reflected;
        int o = configuration;

        for(k=0; k<4; ++k)
        {
            if(cellConfigurations[configuration] != cellConfigurations[o])
            {
                cerr << configuration << ' ' << o << endl;
            }
            if(cellConfigurations[configuration] != cellConfigurations[r])
            {
                cerr << configuration << ' ' << r << endl;
            }

            int i2 = ((o << 2) | (o >> 6)) & 0xff;
            o = i2;

            i2 = ((r << 2) | (r >> 6)) & 0xff;
            r = i2;
        }
    }
}
