#include "chi_runtime.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include "UnColFAlgo/uncolfalgo.h"

int main(int argc, char* argv[])
{
    ChiTech::Initialize(argc,argv);

    chi_log.Log(LOG_0) << "***** UncolFAlgo *****";
    ChiTech::RunBatch(argc,argv);

    chi_log.Log(LOG_0) << "***** Executing custom code *****";
    UncolFAlgo::Solver uncolfalgo(2,1);
    uncolfalgo.Initialize();
    uncolfalgo.Execute();

    ChiTech::Finalize();
    return 0;
}
