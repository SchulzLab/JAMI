package org.mpii.jami.cmi;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;

/**
 * Created by fuksova on 11/9/15.
 * One instance of this class represents one cell of grid in CMI computation in IterativePartitioning
 * Points are represented not by their original double values but as the order the original value in every dimension
 * Important: All methods are designed only for dimension=3
 */
public class Cube {
    private int [][] coordinates; //coordinates of cube in this form {{x-min,x-max},{y-min,y-max},{z-min,z-max}}
    private int dimension;
    private ArrayList<Integer> points;
    private int depth;  //is used only if maxDeep in IterativePartitioning is set to some limiting value
    private int[] dimensions = new int[]{0,1,2};

    public HashSet<Integer> getPointsAfterZeroSplit() {
        return pointsAfterZeroSplit;
    }

    public void setPointsAfterZeroSplit(HashSet<Integer> pointsAfterZeroSplit) {
        this.pointsAfterZeroSplit = pointsAfterZeroSplit;
    }

    private HashSet<Integer> pointsAfterZeroSplit;

    /**
     *
     * @param coordinates enter coordinates in the order x1-min, x1-max, x2-min, x2-max,..
     * @param depth is used only if maxDeep in IterativePartitioning is set to some limiting value
     */
    public Cube(int [] coordinates,int depth){
        this.dimension=coordinates.length/2;
        this.coordinates=new int[dimension][2];
        for (int i = 0; i < dimension; i++) {
            this.coordinates[i][0]=coordinates[2*i];
            this.coordinates[i][1]=coordinates[2*i+1];

        }
        points=new ArrayList<>();
        this.depth =depth;
    }

    public Cube(int[][] coordinates, int depth, int[] dimensions){
        this.dimension = coordinates.length;
        this.coordinates = coordinates;
        points=new ArrayList<>();
        this.depth = depth;
        this.dimensions = dimensions;
    }

    /**
     * Based on cube bounds (coordinates), computes its center
     * @return
     */
    public int[] findCenter(){
        int []center=new int[dimension];
        for (int i = 0; i < dimension; i++) {
            center[i]=(coordinates[i][0]+coordinates[i][1])/2;
        }
        return center;

    }

    /**
     * Add point to the point list of the cube
     * @param pointOrder points are represented by its index in the original point structure data or inverseSortedIndices
     *                   in IterativePartitioning
     */
    public void addPoint(int pointOrder){
        points.add(pointOrder);
    }

    /**
     * Adds Collection of points into this cube
     * @param newPoints points are represented by its index in the original point structure data or inverseSortedIndices
     *                   in IterativePartitioning
     */
    public void addPoints(Collection<Integer> newPoints){
        points.addAll(newPoints);
    }

    /**
     *
     * @return list of points in cube
     */
    public ArrayList<Integer> getPoints() {
        return points;
    }

    /**
     *
     * @return cube coordinates (bounds)
     */
    public int[][] getCoordinates() {
        return coordinates;
    }

    public int[][] cloneCoordinates() {
        int length = coordinates.length;
        int[][] target = new int[length][coordinates[0].length];
        for (int i = 0; i < length; i++) {
            System.arraycopy(coordinates[i], 0, target[i], 0, coordinates[i].length);
        }
        return target;
    }

    public int getDepth() {
        return depth;
    }

    public int getDimension() {
        return dimension;
    }

    public int[] getDimensions() {
        return dimensions;
    }

    public void setDimensions(int[] dimensions) {
        this.dimensions = dimensions;
    }

    public void setCoordinates(int[][] coordinates) {
        this.coordinates = coordinates;
    }
}
