#include <iostream.h>
#include <string>
#include "configurations.h"

int main()
{
    int i = 0;
    
    cerr << "Test: " << cellConfigurations[0] << ' ' << 
                        cellConfigurations[255] << endl;
    
    for(; i<256; ++i)
    {
        int k;
        
        int bits[8];
        for(k=0; k<8; ++k)
        {
            bits[k] = (i >> k) & 1;
        }
        
        int ir  = bits[0] | 
                 (bits[7] << 1) |
                 (bits[6] << 2) |
                 (bits[5] << 3) |
                 (bits[4] << 4) |
                 (bits[3] << 5) |
                 (bits[2] << 6) |
                 (bits[1] << 7);
        
        int r = ir;
        int o = i;
        
        for(k=0; k<4; ++k)
        {
            if(cellConfigurations[i] != cellConfigurations[o])
            {
                cerr << i << ' ' << o << endl;
            }
            if(cellConfigurations[i] != cellConfigurations[r])
            {
                cerr << i << ' ' << r << endl;
            }

            int i2 = ((o << 2) | (o >> 6)) & 0xff;
            o = i2;
            
            i2 = ((r << 2) | (r >> 6)) & 0xff;
            r = i2;
        }
    }
}
