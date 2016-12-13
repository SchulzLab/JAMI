package org.mpii.jami;

import java.util.ArrayList;
import java.util.Collection;

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
        for (Integer newPoint : newPoints) {
            points.add(newPoint);
        }
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


    public int getDepth() {
        return depth;
    }
}
