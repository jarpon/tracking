#include <iostream>
#include <componentlabelling.h>
#include "otsuthresholding.h"
//#include "regionanalysis2.h"
//#include "regionanalysis.h"
#include <regionanalysis3d.h>
#include "thresholding.h"
#include "voxelmatrix.h"
#include <volumehistogramexpansion.h>
#include <voxelmatrixdilatation.h>
#include <voxelmatrixerosion.h>
#include <holesfilling.h>
#include <medianfilter.h>
#include <gaussiangradient.h>
#include <sobelgradient.h>
#include <watershedtransform.h>

#include <cmath>


#define TRACE
#include <trace.h>

VoxelMatrix <float> findParticles(const VoxelMatrix<float>& originalVoxelMatrix, VoxelMatrix<float>& nucleusMask,
                            const string& filename, const string& intermediateProcessesDir )
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
  //regionAnalysis.setLabelMatrix( nucleusMaskCopy );
  regionAnalysis.setValueMatrix( copyVoxelMatrix );
  regionAnalysis.run();
  rangeMask.setSize( copyVoxelMatrix.getSize() );
  rangeMask.setZeros();
  regionAnalysis.setOutputMatrix( rangeMask );

  regionAnalysis.outputFillRegions( REGION_FEATURE_CONTRAST );
  //regionAnalysis.outputFillRegions( REGION_FEATURE_MAXIMUM_VALUE );

  //rangeMask.save( intermediateProcessesDir + filename + "-contrast.vm", true );
  rangeMask.save( intermediateProcessesDir + filename + "-contractness.vm", true );

  OtsuThresholding<float> otsuThresholding;

  Vector <float> featureValues = regionAnalysis.computeRegionFeature( REGION_FEATURE_CONTRAST );
  //Vector <float> featureValues = regionAnalysis.computeRegionFeature( REGION_FEATURE_MAXIMUM_VALUE );
  featureValues.sort();
  Vector <unsigned int> histogram = featureValues.histogram( featureValues.min(), 1, floor(featureValues.max())+1 );

  //float threshold = otsuThresholding.computeThreshold( histogram );
  float threshold = otsuThresholding.computeThreshold( rangeMask );
  EVAL(featureValues);
  EVAL(histogram);
  EVAL(threshold);
  regionAnalysis.thresholdRegions( featureValues, threshold );
  regionAnalysis.run();
  int num = regionAnalysis.condenseRegionLabels();

  VoxelMatrix<float> ccsMask2 = regionAnalysis.getLabelMatrix();
  ccsMask2.save( "/home/jarpon/data/" + filename + "1b.vm", true  );
  regionAnalysis.run();

  rangeMask.save( "/home/jarpon/data/" + filename + "2.vm", true  );
  EVAL(num);
  VoxelMatrix<float> ccsMask = regionAnalysis.getLabelMatrix();
  ccsMask.save( "/home/jarpon/data/" + filename + "3.vm", true  );

  Thresholding<float> altThresholding;
  altThresholding.setBackground( 0.0 );
  altThresholding.setThreshold( threshold );
  altThresholding.apply( rangeMask );

  //labeling the image
  ComponentLabelling<float> componentLabelling;
  componentLabelling.apply( rangeMask );

  //to obtain a better filling, it's applied to each 2D slice instead of the complete 3D stack
  HolesFillingf holesFilling;
  //int sizeZ = ccsMask.getSize3();
  for (int k = 0; k < sizeZ; ++k)  holesFilling.apply( ccsMask[k] );


  ccsMask = rangeMask;

//  VoxelMatrix<float> ccsFillMask = ccsMask;
//  ccsMask.fillIt(ccsFillMask);


//  VoxelMatrix<float> structElement;
//  structElement.setSize(3,3,3);
//  structElement.setOnes();

//  VoxelMatrixDilatation<float> voxelDilatation;
//  voxelDilatation.setStructElt( structElement );
//  voxelDilatation.apply( ccsMask );

// //  VoxelMatrix<float> structElement2;
// //  structElement2.setSize(1,1,1);
// //  structElement2.setOnes();

//  VoxelMatrixErosion<float> voxelErosion;
//  voxelErosion.setStructElt( structElement );
//  voxelErosion.apply( ccsMask );

//  for (int k = 0; k < sizeZ; ++k)  holesFilling.apply( ccsMask[k] );

  ccsMask.setVoxelCalibration( originalVoxelMatrix.getVoxelCalibration() );

  //ccsMask.save( chromocentersDir + filename + ".vm", true );

  return ccsMask;
}

