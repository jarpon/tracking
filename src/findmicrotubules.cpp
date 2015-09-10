
#include <componentlabelling.h>
#include <holesfilling.h>
#include "otsuthresholding.h"
#include "thresholding.h"
#include <voxelmatrix.h>
#include <volumehistogramexpansion.h>
//#include <maxfilter.h>
//#include <minfilter.h>
#include <voxelmatrixerosion.h>
#include <voxelmatrixdilatation.h>
#include <image.h>
#include <medianfilter.h>
#include <gaussiangradient.h>
#include <sstream>

#define TRACE
#include <trace.h>
using namespace std;

VoxelMatrix<float> findMicrotubules(const VoxelMatrix<float>& originalVoxelMatrix,
                                  const string& filename, const string& nucleiDir)
{
  VoxelMatrix<float> microtubulesMask = originalVoxelMatrix;

  MedianFilter<float> medianFilter2;
  medianFilter2.setHalfSize( 1.8 );
  medianFilter2.apply( microtubulesMask );


  OtsuThresholding<float> otsuThresholding;
  otsuThresholding.setForegroundIsAbove( true );
  otsuThresholding.setForeground( 1.0 );
  otsuThresholding.setBackground( 0.0 );
  otsuThresholding.apply( microtubulesMask );

  const int size1 = microtubulesMask.getSize1();
  const int size2 = microtubulesMask.getSize2();
  const int size3 = microtubulesMask.getSize3();

  HolesFillingf holesFilling;

  //to obtain a better filling, it's applied to each 2D slice instead of the complete 3D stack
//  for (int k = 0; k < size3; ++k)  holesFilling.apply( microtubulesMask[k] );
  holesFilling.apply( microtubulesMask );

//this was uncommented 12 june
//  VoxelMatrix<float> nucleusFill = microtubulesMask;
//  microtubulesMask.fillIt(nucleusFill);

//  /*process to improve nuclei segmentation when the nucleoli touch the envelope
//   */
//  VoxelMatrix<float> structElement;
//  structElement.setSize(5,5,5);
//  structElement.setOnes();
//  //different structure element due to problems of segmentation (a kind of smoothing)
//  VoxelMatrix<float> structElement2;
//  structElement2.setSize(3,3,3);
//  structElement2.setOnes();

//  VoxelMatrixDilatation<float> voxelDilatation;
//  voxelDilatation.setStructElt( structElement );
//  voxelDilatation.apply( microtubulesMask );
//  VoxelMatrixErosion<float> voxelErosion;
//  voxelErosion.setStructElt( structElement2 );
//  voxelErosion.apply( microtubulesMask );

//  holesFilling.apply( microtubulesMask );
//  microtubulesMask.fillIt(microtubulesMask);

  /* label each region
   */
  ComponentLabelling<float> componentLabelling;
  componentLabelling.apply( microtubulesMask );

  /* check if some label touches the boundary
   * (that means it's truncated)
   */


//  int labelTruncated = 0;

//  for ( int i = 0; i < size1; ++i )
//    for ( int j = 0; j < size2; ++j )
////      for ( int k = 0; k < size3; ++k )
//          if ( i==0 || i==size1-1 || j==0 || j==size2-1 )
//            if ( microtubulesMask(i,j) != 0 )
//            {
//                labelTruncated = microtubulesMask(i,j);
//                EVAL(labelTruncated);
//                for (int ii = 0; ii < size1; ++ii)
//                for (int jj = 0; jj < size2; ++jj)
////                  for (int kk = 0; kk < size3; ++kk)
//                     if ( microtubulesMask(ii,jj) == labelTruncated )
//                         microtubulesMask(ii,jj) = 0;
//            }

  componentLabelling.apply( microtubulesMask );
  const int numComponents = componentLabelling.getNumLabels();

  if ( numComponents == 0 )
  {
    microtubulesMask = originalVoxelMatrix;
    otsuThresholding.setThreshold( otsuThresholding.getThreshold()*0.8 );
    otsuThresholding.apply( microtubulesMask );
  }

  /* check if there are more than one nucleus in the image
   * regarding to the volumes of the different labels
   */

  int numNuclei = 1;

  if ( numComponents > 1 )
  {
    const Vector<unsigned int> histogram = microtubulesMask.histogram( 0, 1, numComponents+1 );
    const Vector<unsigned int> histogram2 = histogram.copy( 1, histogram.getSize()-1 );
    EVAL(histogram2);

    Thresholding<float> thresholding;
    thresholding.setForegroundIsAbove( true );
    thresholding.setBackground( 0.0 );
    thresholding.setForeground( 1.0 );

    microtubulesMask.setVoxelCalibration( originalVoxelMatrix.getVoxelCalibration() );

    VoxelMatrix<float> tempmicrotubulesMask = microtubulesMask;
    for ( int count = 0; count < numComponents; ++count )
    {
        if ( histogram2.max()/histogram2[count] < 4 && histogram2.maxPos() != count && histogram2[count]>5000 )
      //Z  if ( histogram2.max()/histogram2[count] < 4  )
        {
            ostringstream iss;
            iss << numNuclei;
            for (int k = 0; k < size3; ++k)  thresholding.levelSetMask( tempmicrotubulesMask[k], count+1 );

  //this was uncommented 12 june
        //    VoxelMatrix<float> nucleusFill = tempmicrotubulesMask;
        //    tempmicrotubulesMask.fillIt(nucleusFill);
            //microtubulesMask = nucleusFill;


            //process to improve nuclei segmentation when the nucleoli touch the envelope
            VoxelMatrix<float> structElement1, structElement2, structElement3, structElement4, structElement5;
            structElement1.setSize(3,3,9);
            structElement1.setOnes();
            structElement2.setSize(3,3,7);
            structElement2.setOnes();
            structElement3.setSize(9,9,1);
            structElement3.setOnes();
            structElement4.setSize(5,5,1);
            structElement4.setOnes();
            structElement5.setSize(5,5,3);
            structElement5.setOnes();

//            VoxelMatrixErosion<float> voxelErosion1;
//            voxelErosion1.setStructElt( structElement3 );
//            voxelErosion1.apply( microtubulesMask );

//            VoxelMatrixDilatation<float> voxelDilatation1;
//            voxelDilatation1.setStructElt( structElement1 );
//            voxelDilatation1.apply( microtubulesMask );

//            VoxelMatrixErosion<float> voxelErosion2;
//            voxelErosion2.setStructElt( structElement4 );
//            voxelErosion2.apply( microtubulesMask );

//            VoxelMatrixDilatation<float> voxelDilatation2;
//            voxelDilatation2.setStructElt( structElement2 );
//            voxelDilatation2.apply( microtubulesMask );

            for (int k = 0; k < size3; ++k)  holesFilling.apply( tempmicrotubulesMask[k] );

            tempmicrotubulesMask.save ( nucleiDir + filename + "-" + iss.str() + ".vm", true );
            tempmicrotubulesMask = microtubulesMask;
            ++numNuclei;

        }
    }

    //for (int k = 0; k < size3; ++k)
        thresholding.levelSetMask( microtubulesMask, histogram2.maxPos()+1 );

    //process to improve nuclei segmentation when the nucleoli touch the envelope
    VoxelMatrix<float> structElement1, structElement2, structElement3, structElement4, structElement5;
    structElement1.setSize(3,3,1);
    structElement1.setOnes();
    structElement2.setSize(3,3,1);
    structElement2.setOnes();
    structElement3.setSize(9,9,1);
    structElement3.setOnes();
    structElement4.setSize(5,5,1);
    structElement4.setOnes();
    structElement5.setSize(5,5,1);
    structElement5.setOnes();

//    VoxelMatrixErosion<float> voxelErosion1;
//    voxelErosion1.setStructElt( structElement3 );
//    voxelErosion1.apply( microtubulesMask );

//    VoxelMatrixDilatation<float> voxelDilatation1;
//    voxelDilatation1.setStructElt( structElement1 );
//    voxelDilatation1.apply( microtubulesMask );

//    VoxelMatrixErosion<float> voxelErosion2;
//    voxelErosion2.setStructElt( structElement4 );
//    voxelErosion2.apply( microtubulesMask );

//    VoxelMatrixDilatation<float> voxelDilatation2;
//    voxelDilatation2.setStructElt( structElement2 );
//    voxelDilatation2.apply( microtubulesMask );

    for (int k = 0; k < size3; ++k)  holesFilling.apply( microtubulesMask[k] );

    if ( numNuclei != 1 )
    {
        ostringstream iss;
        iss << numNuclei;
        microtubulesMask.save ( nucleiDir + filename + "-" + iss.str() + ".vm", true );
    }
    else
        microtubulesMask.save ( nucleiDir + filename + ".vm", true );
  }

  else if ( numComponents == 1 )
  {
      //process to improve nuclei segmentation when the nucleoli touch the envelope
    VoxelMatrix<float> structElement1, structElement2, structElement3, structElement4, structElement5;
    structElement1.setSize(3,3,1);
    structElement1.setOnes();
    structElement2.setSize(3,3,1);
    structElement2.setOnes();
    structElement3.setSize(9,9,1);
    structElement3.setOnes();
    structElement4.setSize(5,5,1);
    structElement4.setOnes();
    structElement5.setSize(5,5,1);
    structElement5.setOnes();

//    VoxelMatrixErosion<float> voxelErosion1;
//    voxelErosion1.setStructElt( structElement3 );
//    voxelErosion1.apply( microtubulesMask );

//    VoxelMatrixDilatation<float> voxelDilatation1;
//    voxelDilatation1.setStructElt( structElement1 );
//    voxelDilatation1.apply( microtubulesMask );

//    VoxelMatrixErosion<float> voxelErosion2;
//    voxelErosion2.setStructElt( structElement4 );
//    voxelErosion2.apply( microtubulesMask );

//    VoxelMatrixDilatation<float> voxelDilatation2;
//    voxelDilatation2.setStructElt( structElement2 );
//    voxelDilatation2.apply( microtubulesMask );

      //for (int k = 0; k < size3; ++k)
        holesFilling.apply( microtubulesMask);

      microtubulesMask.save ( nucleiDir + filename + ".vm", true );
  }

  //another 3D filling
  //for (int k = 0; k < size3; ++k)  holesFilling.apply( microtubulesMask[k] );


  return microtubulesMask;
}

