#include "dsmc_io.h"

DSMC_IO::DSMC_IO(Settings *settings_, System *system_) {
    settings = settings_;
    system = system_;
}

void DSMC_IO::state_to_file_binary() {
    int N = system->N;

}
