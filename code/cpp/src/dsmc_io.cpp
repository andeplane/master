#include <dsmc_io.h>
#include <fstream>

using namespace std;

DSMC_IO::DSMC_IO(Settings *settings_, System *system_) {
    settings = settings_;
    system = system_;
}

void DSMC_IO::state_to_file_binary() {
    int N = system->N;
    ofstream myFile ("state", ios::out | ios::binary);
    myFile.write (reinterpret_cast<char*>(N), sizeof(int));

    myFile.write (reinterpret_cast<char*>(N), sizeof(int));

}
