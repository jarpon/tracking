//#include <componentlabelling.h>
#include <dataset.h>
//#include "regionanalysis2.h"
//#include "regionanalysis.h"
#include <cmath>
//#include <marchingcubes.h>
//#include <thresholding.h>
#include <trimesh.h>
#include <sstream>
#include <fileinfo.h>

#define TRACE
#include <trace.h>

void particlesInterdistances(const string& filename,
                           const string& analysisDir, DataSet& individualChromocentersDataset)
{
  const DataSet datasetMicrotubule( analysisDir + filename + "_chromocenters.csv" );
  const int numCompartments = datasetMicrotubule.size()[0];

  Vertices<float> vertices ( 3, numCompartments, 0, 0 );
  for ( int i = 0; i < numCompartments; ++i )
  {
    vertices[i][0] = datasetMicrotubule.getValue<float>( "centroidCoordX", i );
    vertices[i][1] = datasetMicrotubule.getValue<float>( "centroidCoordY", i );
    vertices[i][2] = datasetMicrotubule.getValue<float>( "centroidCoordZ", i );
    EVAL(vertices[i]);
  }

  float tempDistance;
  Vector<float> temp1, temp2;

  for (int i = 0; i < numCompartments; ++i)
  {
    temp1 = vertices[i];

    for (int j = 0; j < numCompartments; ++j)
    {
      if ( i == j )
        tempDistance = 0;
      else
      {
        temp2 = vertices[j];
        tempDistance = temp1.distance(temp2);
      }
      ostringstream iss; //we suppose as much 99 labels
      iss << (j+1);
      individualChromocentersDataset.setValue ( iss.str(), i, tempDistance );
    }
  }
}


void particlesInterdistances(const DataSet& list,
                           const string& analysisDir, DataSet& particlesDistanceDataset)
{

//  particlesDistanceDataset.setValues ( "dataFrame", 0 );
//  particlesDistanceDataset.setValues ( "particle 1 .x", 0 );
//  particlesDistanceDataset.setValues ( "particle 1 .y", 0 );
//  particlesDistanceDataset.setValues ( "particle 1 - movement", 0 );
//  particlesDistanceDataset.setValues ( "particle 1 .x - movement", 0 );
//  particlesDistanceDataset.setValues ( "particle 1 .y - movement", 0 );

  string filename;
  float tempDistance;
  Vector<float> temp1, temp2;

  for ( int k = 0; k < list.size()[0]; ++k )
  {
    Vertices<float> tempVertices ( 2, 1 );

    ostringstream segmented;

    filename = list.getValue<string>( "dataFrame", k );
    EVAL(filename);
    particlesDistanceDataset.setValue ( "dataFrame", k, filename );

    FileInfo fileInfo( analysisDir + filename + "_chromocenters.csv" );


    if ( fileInfo.fileExists() == true )
    {
      segmented << "yes";
      particlesDistanceDataset.setValue ( "segmented", k, segmented.str() );

      const DataSet datasetMicrotubule( analysisDir + filename + "_chromocenters.csv" );
      const int numCompartments = datasetMicrotubule.size()[0];
      EVAL(numCompartments);

      Vertices<float> vertices ( 2, 1 );


      for ( int i = 0; i < numCompartments; ++i )
      {
        vertices[i][0] = datasetMicrotubule.getValue<float>( "centroidCoordX", i );
        vertices[i][1] = datasetMicrotubule.getValue<float>( "centroidCoordY", i );
  //      vertices[i][2] = datasetMicrotubule.getValue<float>( "centroidCoordZ", i );

        EVAL(vertices[i]);

        if (  particlesDistanceDataset.getValue<string>( "segmented", k-1 ) == "yes" )
        {
          ostringstream iss;
          iss << i;

          particlesDistanceDataset.setValue ( "particle " + iss.str() + " .x", k, vertices[i][0] );
          particlesDistanceDataset.setValue ( "particle " + iss.str() + " .y", k, vertices[i][1] );
    //      particlesDistanceDataset.setValue ( "particle " + i + " .z", k, vertices[i][2] );

          temp2 = vertices[i];
          temp1 = tempVertices[i];

          tempDistance = temp2.distance(temp1);

          particlesDistanceDataset.setValue ( "particle " + iss.str() + " - movement", k, tempDistance );

          particlesDistanceDataset.setValue ( "particle " + iss.str() + " .x - movement", k, temp2[0] - temp1[0] );
          particlesDistanceDataset.setValue ( "particle " + iss.str() + " .y - movement", k, temp2[1] - temp1[1] );
    //      particlesDistanceDataset.setValue ( "particle " + i + " .z - movement", k, temp2[2] - temp1[2] );

        }
        else
        {
          ostringstream iss;
          iss << i;

          particlesDistanceDataset.setValue ( "particle " + iss.str() + " .x", k, vertices[i][0] );
          particlesDistanceDataset.setValue ( "particle " + iss.str() + " .y", k, vertices[i][1] );
    //      particlesDistanceDataset.setValue ( "particle " + i + " .z", k, vertices[i][2] );

        }

      }

      tempVertices = vertices;

    }

    else
    {
      segmented << "no";
      particlesDistanceDataset.setValue ( "segmented", k, segmented.str() );
//      particlesDistanceDataset.setValue ( "particle 1 .x", k, 0 );
//      particlesDistanceDataset.setValue ( "particle 1 .y", k, 0 );
////      particlesDistanceDataset.setValue ( "particle " + i + " .z", k, vertices[i][2] );

//      particlesDistanceDataset.setValue ( "particle 1 - movement", k, 0 );

//      particlesDistanceDataset.setValue ( "particle 1 .x - movement", k, 0 );
//      particlesDistanceDataset.setValue ( "particle 1 .y - movement", k, 0 );
////      particlesDistanceDataset.setValue ( "particle " + i + " .z - movement", k, temp2[2] - temp1[2] );


      tempVertices.setSize( 0, true );
    }

  }


}

