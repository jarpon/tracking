#include <iostream>
#include <dataset.h>
#include <dirent.h>
#include <errno.h>
#include <fileinfo.h>
//#include <QString>
#include <unistd.h>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <voxelmatrix.h>
#include <trimesh.h>
#include <sstream>
#include <programerror.h>
#include <stopwatch.h>

#define TRACE
#include <trace.h>
using namespace std;

extern VoxelMatrix<float> findMicrotubules(const VoxelMatrix<float>&,
                                         const string&, const string&);
extern VoxelMatrix<float> findParticles( const VoxelMatrix<float>&, VoxelMatrix<float>&,
                                   const string&, const string&);
extern VoxelMatrix<float> findParticlesManually( const VoxelMatrix<float>&, VoxelMatrix<float>&,
                                   const string&, const string&, const string&);
extern void particlesAnalysis(VoxelMatrix<float>&, const string&, const string&,
                                  const int&, int&, DataSet&, DataSet&, DataSet&);
extern void particlesInterdistances(const DataSet&,
                           const string&, DataSet&);

int main(int argc, char* argv[])
{

  string filepath, filename, parentDir, originalDir, originalVMDir, nucleiDir, chromocentersDir, intermediateProcessesDir, analysisDir, shapesDir;

  if ( (argc == 1) || ( argv[1] == std::string("-h") ) || ( argv[1] == std::string("--help") ) )
  {
    // Check the value of argc. If not enough parameters have been passed, inform user and exit.
    // Inform the user of how to use the program
    cout << "                                                         " << endl;
    cout << "Usage: tracking [OPTION process to apply] FILES... " << endl;
  }

  else if ( argv[1] == std::string("-p") &&
     ( argv[2] == std::string("1") || argv[2] == std::string("3") || argv[2] == std::string("3m") ||
       argv[2] == std::string("4") || argv[2] == std::string("5") ) && ( argc > 3 ) )
  {
    Stopwatch stopWatch;
    stopWatch.start( "Global process" );

    const string process = argv[2];
    EVAL(process);

    DataSet nucleiDataset;
    DataSet particlesDataset;

    //these are just counters to use at the datafiles
    int numNucleus = 0;
    int totalNumParticles = 0;

    DataSet framesList;

    for (int i = 3; i < argc; i++)
    {
      filename = argv[i];

      FileInfo fileInfo (filename);

      if ( fileInfo.isAbsolutePath() == false )
      {
        originalDir = getcwd( argv[i], 4048 );
        originalDir = originalDir + "/";
        filepath = originalDir + filename;
        FileInfo fileInfo (filepath);
        parentDir = originalDir.substr(0,originalDir.find_last_of("/\\"));
        parentDir = parentDir.substr(0,parentDir.find_last_of("/\\")) + "/";
        filename = fileInfo.baseName();
      }
      else
      {
        filepath = filename;
        filename = fileInfo.baseName();
      }

      framesList.setValue( "dataFrame", i-3, filename );
     // EVAL(filename);
    }

    for (int i = 3; i < argc; i++)
    {
      filename = argv[i];

      // only vm files are processed
      if(filename.substr(filename.find_last_of(".") + 1) == "vm")
      {
//        DataSet individualChromocentersDataset;


        EVAL(filename);
        FileInfo fileInfo (filename);

        if ( fileInfo.isAbsolutePath() == false )
        {
          originalDir = getcwd( argv[i], 4048 );
          originalDir = originalDir + "/";
          filepath = originalDir + filename;
          FileInfo fileInfo (filepath);
          parentDir = originalDir.substr(0,originalDir.find_last_of("/\\"));
          parentDir = parentDir.substr(0,parentDir.find_last_of("/\\")) + "/";
          filename = fileInfo.baseName();
        }

        else
        {
          filepath = filename;
          filename = fileInfo.baseName();
          originalDir = fileInfo.dirName();
          parentDir = originalDir.substr(0,originalDir.find_last_of("/\\"));
          parentDir = parentDir.substr(0,parentDir.find_last_of("/\\"));
          //originalDir = parentDir + "/originals_vm/";
        }

        nucleiDir = parentDir + "/segmented_nuclei/";
        originalVMDir = parentDir + "/originals_vm/";
        intermediateProcessesDir = parentDir + "/intermediate_processes/";
        chromocentersDir = parentDir + "/segmented_chromocenters/";
        analysisDir = parentDir + "/analysis/";

        EVAL(filename);
        EVAL(originalVMDir);

        if ( process == "1" )
        {
          ENTER("Microtubules segmentation, looking for more than 1");
          VoxelMatrix<float> originalVoxelMatrix( originalVMDir + filename + ".vm" );
          VoxelMatrix<float> nucleusMask;
//          PixelMatrix<float> originalPixelMatrix;
//          originalPixelMatrix.convertFromImage( originalVoxelMatrix );
          nucleusMask = findMicrotubules( originalVoxelMatrix, filename, nucleiDir );
          //nucleusMask.save ( nucleiDir + filename + ".vm", true );
          LEAVE();
        }

        else if ( process == "3" )
        {
          ENTER("Chromocenters segmentation");
          string originalName = filename.substr( 0,filename.find_last_of("-")  );
          VoxelMatrix<float> originalVoxelMatrix( originalVMDir + originalName + ".vm" );
          //VoxelMatrix<float> originalVoxelMatrix( originalVMDir + filename + ".vm" );
          VoxelMatrix<float> nucleusMask ( nucleiDir + filename + ".vm" );
          VoxelMatrix<float> ccsMask;
          ccsMask = findParticles( originalVoxelMatrix, nucleusMask, filename, intermediateProcessesDir);
          ccsMask.save ( chromocentersDir + filename + ".vm", true );
          LEAVE();
        }

        else if ( process == "3m" )
        {
          ENTER("Chromocenters manual segmentation");
          string originalName = filename.substr( 0,filename.find_last_of("-")  );
          VoxelMatrix<float> originalVoxelMatrix( originalVMDir + originalName + ".vm" );
          //VoxelMatrix<float> originalVoxelMatrix( originalVMDir + filename + ".vm" );
          VoxelMatrix<float> nucleusMask ( nucleiDir + filename + ".vm" );
          VoxelMatrix<float> ccsMask;
          ccsMask = findParticlesManually( originalVoxelMatrix, nucleusMask, filename, intermediateProcessesDir, chromocentersDir);
          ccsMask.save ( chromocentersDir + filename + ".vm", true );
          LEAVE();
        }
        
        else if ( process == "4" )
        {
          ENTER("Particles quantification");
          DataSet individualParticlesDataset;
          VoxelMatrix<float> ccsMask ( chromocentersDir + filename + ".vm" );
          particlesAnalysis( ccsMask, filename, parentDir, numNucleus, totalNumParticles, nucleiDataset, individualParticlesDataset, particlesDataset );
          individualParticlesDataset.save(analysisDir + filename + "_chromocenters.csv", true );
          ++numNucleus;
          LEAVE();
        }

        else if ( process == "5" )
        {
          ENTER("Particles tracking");
          DataSet particlesDistanceDataset;
          particlesInterdistances( framesList, analysisDir, particlesDistanceDataset );
          particlesDistanceDataset.save(analysisDir + "particlesTracking.csv", true );
          LEAVE();
        }
      }

      else cout << "Error opening the image, it must be a VoxelMatrix image " << endl;

      if ( process == "5" )
        return 0;

    }

    if ( argv[2] == std::string("4") )
    {
      nucleiDataset.save(analysisDir + "nuclei_extended.csv", true );
      particlesDataset.save(analysisDir + "particles.csv", true );
    }

    stopWatch.stop( "Global process" );
    stopWatch.print();
  }

  return 0;
}
