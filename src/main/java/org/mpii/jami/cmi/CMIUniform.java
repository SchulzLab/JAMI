package org.mpii.jami.cmi;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by fuksova on 1/21/16.
 * Compute conditional mutual information on uniform grid
 */
public class CMIUniform {
    ArrayList<List<Double>> data;  //Input data, double[] arrays have to have all the same length
    int pointsSize;   //length of double[] arrays in data
    int dimension;    //size of data. Currently the dimension has to be equal to 3!




    public CMIUniform(ArrayList<List<Double>> data){
        this.data=data;
        pointsSize=data.get(0).size();
        dimension=data.size();

    }

    /**
     * Returns value of CMI for data I(data.get(0),data.get(1)|data.get(2))
     * @param numberOfBins Gives number of bins in each dimension for uniform the grid where we compute CMI,
     * @return value of CMI
     */
    public double computeCMI(int numberOfBins){

        double [][]minMax=new double[dimension][2];

        for (int i = 0; i < dimension; i++) {
            minMax[i]=findMinAndMax(data.get(i));
        }
        double[]binWidth=new double[dimension];
        for (int i = 0; i < dimension; i++) {
            binWidth[i]=(minMax[i][1]-minMax[i][0])/numberOfBins;
            if(binWidth[i]==0){
                return 0;
            }
        }

        CubeDoubleBounds[][][] allCubes=new CubeDoubleBounds[numberOfBins][numberOfBins][numberOfBins];
        double[]currentCoord=new double[dimension];
        for (int i = 0; i < numberOfBins; i++) {
            currentCoord[0]=i;
            for (int j = 0; j < numberOfBins; j++) {
                currentCoord[1]=j;
                for (int k = 0; k < numberOfBins; k++) {
                    currentCoord[2]=k;
                    double[][]coordinates=new double[dimension][2];
                    for (int l = 0; l < dimension; l++) {
                        coordinates[l][0]=currentCoord[l]*binWidth[l]+minMax[l][0];
                        if(currentCoord[l]==numberOfBins-1){
                            coordinates[l][1]=minMax[l][0]+(currentCoord[l]+1)*binWidth[l]+binWidth[l]*0.01;
                        }
                        else{
                            coordinates[l][1]=minMax[l][0]+(currentCoord[l]+1)*binWidth[l];
                        }
                    }
                    allCubes[i][j][k]=new CubeDoubleBounds(coordinates);
                }
            }

        }
        for (int i = 0; i < pointsSize; i++) {
            int[] cubeOrder=new int[dimension];
            for (int j = 0; j < dimension; j++) {
                double coord=data.get(j).get(i);
                cubeOrder[j]=(int)Math.min(Math.floor((coord-minMax[j][0])/binWidth[j]),numberOfBins-1);
            }
            CubeDoubleBounds selectedCube = allCubes[cubeOrder[0]][cubeOrder[1]][cubeOrder[2]];
            selectedCube.addPoint(i);
            //for debugging purposes:
            for (int j = 0; j < dimension; j++) {
                double coord=data.get(j).get(i);
                if(coord<selectedCube.getCoordinates()[j][0]||coord>=selectedCube.getCoordinates()[j][1]){
                    System.out.println("Wrong point assignment point "+i);
                }
            }
        }

        double cmi=0.0;
        for (int i = 0; i < numberOfBins; i++) {
            for (int j = 0; j < numberOfBins; j++) {
                for (int k = 0; k < numberOfBins; k++) {
                    cmi+=computeCMIInCube(allCubes,i,j,k,numberOfBins);
                }
            }
        }
        return cmi/pointsSize;

    }

    /**
     * Finds maximal and minimal value in array
     * @param array
     * @return returns double[] of length 2, first value is minimum, second is maximum
     */
    private double[] findMinAndMax(List<Double> array){
        double min=Double.MAX_VALUE;
        double max=Double.MIN_VALUE;
        for (int i = 0; i < array.size(); i++) {
            if(array.get(i)<min){
                min=array.get(i);
            }
            if(array.get(i)>max){
                max=array.get(i);
            }
        }
        double[]result=new double[2];
        result[0]=min;
        result[1]=max;
        return result;
    }

    /**
     * Computes CMI in one cube given by allCubes[xOrder][yOrder][zOrder]
     * @param allCubes 3D array of all cubes (cells) in uniform grid
     * @param numberOfBins number of bins in each dimension in uniform grid
     * @return value of CMI
     */
    private double computeCMIInCube(CubeDoubleBounds[][][] allCubes,int xOrder,int yOrder,int zOrder,int numberOfBins){
        int nz=0;
        double cmi=0.0;
        int sizePointsInCube=allCubes[xOrder][yOrder][zOrder].getPoints().size();
        for (int i = 0; i < numberOfBins; i++) {
            for (int j = 0; j < numberOfBins; j++) {
                nz+=allCubes[i][j][zOrder].getPoints().size();
            }
        }
        double nxz=0.0;
        for (int i = 0; i < numberOfBins; i++) {
            nxz+=allCubes[xOrder][i][zOrder].getPoints().size();
        }
        double nyz=0.0;
        for (int i = 0; i < numberOfBins; i++) {
            nyz+=allCubes[i][yOrder][zOrder].getPoints().size();
        }
        double cond=1;
        if(nxz*nyz>0&&nz*sizePointsInCube>0){
            cond=((double)nz*sizePointsInCube)/(nxz*nyz);
        }
        cmi=sizePointsInCube*Math.log(cond);

        return cmi;

    }
}
