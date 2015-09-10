#include <iostream>
#include <componentlabelling.h>
#include <regionanalysis3d.h>
#include "thresholding.h"
#include "voxelmatrix.h"
#include <voxelmatrixdilatation.h>
#include <voxelmatrixerosion.h>
#include <holesfilling.h>
#include <medianfilter.h>
#include <gaussiangradient.h>
#include <sobelgradient.h>
#include <watershedtransform.h>

#include <cmath>
#include <string>
#include <sstream>

using namespace std;

#define TRACE
#include <trace.h>

void findParticlesManually(const VoxelMatrix<float>& originalVoxelMatrix, VoxelMatrix<float>& nucleusMask,
                            const string& filename, const string& intermediateProcessesDir , const string& chromocentersDir)
{
  VoxelMatrix<float> gradientMatrix = originalVoxelMatrix;
  VoxelMatrix<float> nucleusMaskCopy = nucleusMask;

//  MedianFilter<float> medianFilter;
//  medianFilter.setHalfSize( 2 );
//  medianFilter.setNumIterations( 2 );
//  medianFilter.apply( gradientMatrix );

  int sizeZ = originalVoxelMatrix.getSize3();
  if ( nucleusMaskCopy.max().max().max() > 1 )
  {
      Thresholding<float> thresholding;
      thresholding.setBackground( 0.0 );
      thresholding.setForeground( 1.0 );
      thresholding.setThreshold( 1.0 );
      thresholding.apply( nucleusMaskCopy );
  }

  GaussianGradient<float> gaussianGradient;
  gaussianGradient.MaskVoxelMatrixProcessing<float>::setMask( nucleusMaskCopy );
  gaussianGradient.setSigma( 2 );

  for (int k = 0; k < sizeZ; ++k) gaussianGradient.apply( gradientMatrix [k] );
  //gaussianGradient.apply( gradientMatrix );
  gradientMatrix.save( intermediateProcessesDir + filename + "-gradient.vm", true );
  //gradientMatrix.operator /=( 16 );

  VoxelMatrix <float> regionMatrix = gradientMatrix;
  WatershedTransform<float> watershedTransform;
  watershedTransform.MaskVoxelMatrixProcessing<float>::setMask( nucleusMaskCopy );
  watershedTransform.apply( regionMatrix );

  EVAL("1");
//  EVAL(regionMatrix.getSize());

  VoxelMatrix<float> rangeMask, copyVoxelMatrix = originalVoxelMatrix;

  RegionAnalysis3D<float> regionAnalysis;
  regionAnalysis.setLabelMatrix( regionMatrix );
  regionAnalysis.setValueMatrix( copyVoxelMatrix );
  regionAnalysis.run();

  rangeMask.setSize( copyVoxelMatrix.getSize() );
  rangeMask.setZeros();
  regionAnalysis.setOutputMatrix( rangeMask );

  regionAnalysis.outputFillRegions( REGION_FEATURE_CONTRAST );


  Vector <float> featureValues = regionAnalysis.computeRegionFeature( REGION_FEATURE_CONTRAST );


  //const Vector <float> featureValues = regionAnalysis.allRegionValues( rangeMask );

  float inputThreshold;
  string input = "";
  while (true) {
   cout << "Please enter the required value for the correct segmentation:\n>";
   getline(cin, input);

   // This code converts from string to number safely.
   stringstream myStream(input);
   if (myStream >> inputThreshold)
     break;
   cout << "Invalid number, please try again" << endl;
  }
  cout << "You are going to use: " << inputThreshold << " as threshold" << endl << endl;

  EVAL(inputThreshold);

  if ( inputThreshold == -1 )
  {
    cout << "You have discarded this stack" << endl;
    return;
  }
  else
  {

//    Thresholding<float> thresholding;
//    thresholding.setBackground( 0.0 );
//    thresholding.setThreshold( inputThreshold );
//    thresholding.applyAlternative( ccsMask );

  //  HolesFillingf holesFilling;
  //  int sizeZ = originalVoxelMatrix.getSize3();

    //to obtain a better filling, it's applied to each 2D slice instead of the complete 3D stack
  //  for (int k = 0; k < sizeZ; ++k)  holesFilling.apply( ccsMask[k] );
  //  holesFilling.apply( ccsMask );

    regionAnalysis.thresholdRegions( featureValues, inputThreshold );
    regionAnalysis.run();
    int num = regionAnalysis.condenseRegionLabels();
    regionAnalysis.run();

    EVAL(num);
    VoxelMatrix<float> ccsMask = regionAnalysis.getLabelMatrix();

    VoxelMatrix<float> structElement;
    structElement.setSize(3,3,3);
    structElement.setOnes();

    VoxelMatrixDilatation<float> voxelDilatation;
    voxelDilatation.setStructElt( structElement );
    voxelDilatation.apply( ccsMask );

   //  VoxelMatrix<float> structElement2;
   //  structElement2.setSize(1,1,1);
   //  structElement2.setOnes();

    VoxelMatrixErosion<float> voxelErosion;
    voxelErosion.setStructElt( structElement );
    voxelErosion.apply( ccsMask );

    ccsMask.setVoxelCalibration( originalVoxelMatrix.getVoxelCalibration() );

    ccsMask.save ( chromocentersDir + filename + ".vm", true );

  }

  //return ccsMask;
  return;
}

