#include "image.h"

mat Image::readBMP(const char* filename)
{
    int i;
    FILE* f = fopen(filename, "rb");
    unsigned char info[54];
    fread(info, sizeof(unsigned char), 54, f); // read the 54-byte header

    // extract image height and width from header
    int fileSize = *(int*)&info[2];

    int width  = *(int*)&info[18];
    int height = *(int*)&info[22];

    int pixels = width*height;

    int bytesLeft = fileSize-54; // 54 bytes for the header, we need the rest

    unsigned char* data = new unsigned char[bytesLeft];
    fread(data, sizeof(unsigned char), bytesLeft, f); // read the rest of the data at once
    fclose(f);

    mat img = zeros<mat>(width,height);

    int pixelCount = 0;
    int pixelIndex = 0;

    int row = 0; // BMP-files start with the lower left pixel, row for row
    int col = 0;

    for(i = 0; i < pixels; i++)
    {
        int r = data[pixelIndex+2]; // RGB values are interchanged in the file format
        int g = data[pixelIndex+1];
        int b = data[pixelIndex];

        double avg = 1.0*(r+g+b)/3.0/255.0; // If we have black/white only, the average is 0 (black) to 255 (white)
        img(col,height - row - 1) = avg;

        pixelCount++;
        pixelIndex+=3; // Each pixel has 3 bytes, RGB

        if(++col == width) {
            // BMP-format is stupid since it wants each row to have 4n bytes, so it adds
            // the remaining pixels before next row
            int padding = col % 4;

            col = 0;
            pixelIndex+=padding;
            row++;
        }
    }

    return img.t();
}

Image::Image()
{

}

mat Image::loadBMP(char *filename) {
    img = readBMP(filename);
    return img;
}
