#include <componentlabelling.h>
#include <dataset.h>
//#include "regionanalysis2.h"
//#include "regionanalysis.h"
#include <regionanalysis3d.h>
#include <regionanalysis2d.h>
#include <cmath>
#include <marchingcubes.h>
#include <thresholding.h>
#include <trimesh.h>

#define TRACE
#include <trace.h>

void particlesAnalysis(VoxelMatrix<float>& ccsMask, const string& filename, const string& parentDir,
                           const int& numNucleus, int& totalNumCCs,
                           DataSet& nucleiDataset, DataSet& individualParticlesDataset, DataSet& particlesDataset)
{
  string originalName = filename.substr( 0,filename.find_last_of("-")  );

  VoxelMatrix<float> originalVoxelMatrix( parentDir + "/originals_vm/" + originalName + ".vm" );

  RegionAnalysis3D<float> regionAnalysisParticles;
  regionAnalysisParticles.setLabelMatrix( ccsMask );
  regionAnalysisParticles.setValueMatrix( originalVoxelMatrix );
  regionAnalysisParticles.run();
  EVAL (regionAnalysisParticles.numRegions() );

  int num = regionAnalysisParticles.condenseRegionLabels();
  regionAnalysisParticles.run();

  //TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + "_nucleus.tm" );
  TriMesh<float> nucleusTriMesh ( parentDir + "/shapes/" + filename + ".tm" );

//  PixelMatrix<float> slice;
//  slice = originalVoxelMatrix[0];


  RegionAnalysis2D<float> regionAnalysisParticles2D;
  regionAnalysisParticles2D.setLabelMatrix( ccsMask[0] );
  regionAnalysisParticles2D.setValueMatrix( originalVoxelMatrix[0] );
  regionAnalysisParticles2D.run();
  EVAL (regionAnalysisParticles2D.numRegions() );

//  int num2 = regionAnalysisParticles2D.condenseRegionLabels();
//  regionAnalysisParticles2D.run();


  Vector<float> centroid(3);
  Vector<float> vertexTriMesh(3);

  // get the name of the class
  string classif = parentDir;
  classif = classif.substr(classif.find_last_of("/\\")+1,classif.length());

  Vertices<float> centroids = regionAnalysisParticles.regionCentroids();
  Vector<float> ccsVolume = regionAnalysisParticles.computeRegionFeature( REGION_FEATURE_VOLUME );
  Vector<float> ccsEqRadius = regionAnalysisParticles.computeRegionFeature( REGION_FEATURE_EQUIVALENT_RADIUS );
  Vector<float> ccsSurfaceArea = regionAnalysisParticles.computeRegionFeature( REGION_FEATURE_SURFACE_AREA );
  Vector<float> ccsArea = regionAnalysisParticles2D.computeRegionFeature( REGION_FEATURE_AREA );
//  Vector<float> ccsRelativeVolume = ccsVolume.operator *( nucleusVolume ); //relative volume of cc within the nucleus
  Vector<float> ccsFlatness = regionAnalysisParticles.computeRegionFeature( REGION_FEATURE_FLATNESS );
  Vector<float> ccsElongation = regionAnalysisParticles.computeRegionFeature( REGION_FEATURE_ELONGATION );
  Vector<float> ccsSphericity = regionAnalysisParticles2D.computeRegionFeature( REGION_FEATURE_COMPACTNESS );
  Vector<float> ccsIntegratedDensity = regionAnalysisParticles.computeRegionFeature( REGION_FEATURE_INTEGRATED_INTENSITY );

  nucleiDataset.setValue ( "name", numNucleus, filename );//filename
  nucleiDataset.setValue ( "class", numNucleus, classif );//classification: mutant, tissue, etc.
  nucleiDataset.setValue ( "ccsNumber", numNucleus, regionAnalysisParticles.numRegions() );//number of ccs obtained in the nucleus
  nucleiDataset.setValue ( "ccsVolume", numNucleus, ccsVolume.sum() );//total ccs volume
  nucleiDataset.setValue ( "ccsVolume", numNucleus, ccsVolume.sum() );//total ccs volume
  //  nucleiDataset.setValue ( "volRHF", numNucleus, ccsVolume.sum() / nucleusVolume  );//relative volume of ccs regarding the complete nucleus
  nucleiDataset.setValue ( "ccsIntensity", numNucleus, ccsIntegratedDensity.sum()*ccsVolume.sum() );//total absolute intensity of ccs
  nucleiDataset.setValue ( "ccsIntegratedDensity", numNucleus, ccsIntegratedDensity.sum() );//integrated density of ccs taking into account real volume
//  nucleiDataset.setValue ( "intRHF", numNucleus, ( ccsIntegratedDensity.sum()*ccsVolume.sum() ) / nucleusIntensity );//RHF: rate of heterochromatin (int.Density of ccs/ int.Density of nuclei)



  EVAL ( num );

/**///chromocenters individual information
  for (int numCC = 0; numCC < regionAnalysisParticles.numRegions(); numCC++ )
  {
    particlesDataset.setValue ( "name", numCC+totalNumCCs, filename );
    particlesDataset.setValue ( "class", numCC+totalNumCCs, classif );//classification: mutant, tissue, etc.
    particlesDataset.setValue ( "idCC", numCC+totalNumCCs, numCC+1 );

    centroid = centroids[numCC];
    particlesDataset.setValue ( "centroidCoordX", numCC+totalNumCCs, centroid[X] );
    particlesDataset.setValue ( "centroidCoordY", numCC+totalNumCCs, centroid[Y] );
    particlesDataset.setValue ( "centroidCoordZ", numCC+totalNumCCs, centroid[Z] );

    MarchingCubes<float> marchingCubes;
    TriMesh<float> triMesh;
    VoxelMatrix<float> currentLabeledVM = ccsMask;
    Thresholding<float> thresholding;
    thresholding.setForeground( 1.0 );
    thresholding.setBackground( 0.0 );
    thresholding.levelSetMask( currentLabeledVM, numCC+1 );
    triMesh = marchingCubes.buildMesh( currentLabeledVM, 0.5, true );
    triMesh.scale( originalVoxelMatrix.getVoxelCalibration().getVoxelSize() );

    nucleusTriMesh.closestPoint( centroid, vertexTriMesh );
    float distanceToBorder = centroid.distance( vertexTriMesh );

    particlesDataset.setValue ( "equivalentRadius_vm", numCC+totalNumCCs, ccsEqRadius[numCC] );
    particlesDataset.setValue ( "ccVolume_vm", numCC+totalNumCCs, ccsVolume[numCC] );
    particlesDataset.setValue ( "ccArea", numCC+totalNumCCs, ccsArea[numCC] );
    particlesDataset.setValue ( "distanceToTheBorder", numCC+totalNumCCs, distanceToBorder );
    particlesDataset.setValue ( "voxelSizeUnit", numCC+totalNumCCs, originalVoxelMatrix.getVoxelCalibration().getLengthUnit().symbol() + "^3" );//real voxel size unit
//    particlesDataset.setValue ( "ccRelativeVolume", numCC+totalNumCCs, ccsRelativeVolume[numCC] );
    particlesDataset.setValue ( "flatness", numCC+totalNumCCs, ccsFlatness[numCC] );
    particlesDataset.setValue ( "elongation", numCC+totalNumCCs, ccsElongation[numCC] );
    particlesDataset.setValue ( "circularity", numCC+totalNumCCs, ccsSphericity[numCC] );
    particlesDataset.setValue ( "surfaceArea", numCC+totalNumCCs, ccsSurfaceArea[numCC] );
    particlesDataset.setValue ( "ccsIntensity", numCC+totalNumCCs, ccsIntegratedDensity[numCC] * ccsVolume[numCC] );
    particlesDataset.setValue ( "ccsIntegratedDensity", numCC+totalNumCCs, ccsIntegratedDensity[numCC] );
    //this last relativeCCsIntensity is the intensity of each cc divided by the total cc's intensity
    particlesDataset.setValue ( "relativeCCsIntensity", numCC+totalNumCCs, ccsIntegratedDensity[numCC] * ccsVolume[numCC] / ( ccsIntegratedDensity.sum() * ccsVolume.sum() ) );

    individualParticlesDataset.setValue ( "id", numCC, filename );
    individualParticlesDataset.setValue ( "idCC", numCC, numCC+1 );
    individualParticlesDataset.setValue ( "centroidCoordX", numCC, centroid[X] );
    individualParticlesDataset.setValue ( "centroidCoordY", numCC, centroid[Y] );
    individualParticlesDataset.setValue ( "centroidCoordZ", numCC, centroid[Z] );
    individualParticlesDataset.setValue ( "equivalentRadius_vm", numCC, ccsEqRadius[numCC] );
    individualParticlesDataset.setValue ( "distanceToTheBorder", numCC, distanceToBorder );
    individualParticlesDataset.setValue ( "ccArea", numCC, ccsArea[numCC] );
    individualParticlesDataset.setValue ( "ccVolume_vm", numCC, ccsVolume[numCC] );
    individualParticlesDataset.setValue ( "voxelSizeUnit", numCC, originalVoxelMatrix.getVoxelCalibration().getLengthUnit().symbol() + "^3" );//real voxel size unit
//    individualParticlesDataset.setValue ( "ccRelativeVolume", numCC, ccsRelativeVolume[numCC] );
    individualParticlesDataset.setValue ( "flatness", numCC, ccsFlatness[numCC] );
    individualParticlesDataset.setValue ( "elongation", numCC, ccsElongation[numCC] );
    individualParticlesDataset.setValue ( "circularity", numCC, ccsSphericity[numCC] );
    individualParticlesDataset.setValue ( "surfaceArea", numCC, ccsSurfaceArea[numCC] );
    individualParticlesDataset.setValue ( "ccsIntensity", numCC, ccsIntegratedDensity[numCC] * ccsVolume[numCC] );
    individualParticlesDataset.setValue ( "ccsIntegratedDensity", numCC, ccsIntegratedDensity[numCC] );
    individualParticlesDataset.setValue ( "relativeCCsIntensity", numCC, ccsIntegratedDensity[numCC] * ccsVolume[numCC] / ( ccsIntegratedDensity.sum() * ccsVolume.sum() ) );

  }

  totalNumCCs += regionAnalysisParticles.numRegions();

  //return particlesDataset;
}

