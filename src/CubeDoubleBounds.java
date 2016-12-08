import java.util.ArrayList;
import java.util.Collection;

/**
 * Created by fuksova on 1/21/16.
 * This class is used for computation of conditional mutual information in uniform grid. One CubeDoubleBounds
 * states for one cell in uniform grid
 */
public class CubeDoubleBounds {
    private double [][] coordinates; //bounds of this cube. That is coordinates[0][1], first index means order of coordinate
                                        //second index: 0: lower bound, 1: upper bound
    private int dimension;
    private ArrayList<Integer> points;


    public CubeDoubleBounds(double[][] coordinates){
        this.coordinates=coordinates;
        this.dimension=coordinates.length;
        points=new ArrayList<>();
    }

    /**
     * Adds one point into the list of points in this cube
     * @param point index of point in CMIUniform.data
     */
    public void addPoint(int point){
        points.add(point);
    }

    /**
     * Adds more points into the list of points in this cube
     * @param newPoints Collection of points, one point is represented as its index in CMIUniform.data
     */
    public void addPoints(Collection<Integer> newPoints){
        this.points.addAll(newPoints);
    }

    /**
     *
     * @return list of points in this cube
     */
    public ArrayList<Integer> getPoints() {
        return points;
    }

    /**
     *
     * @return bounds of this cube
     */
    public double[][] getCoordinates() {
        return coordinates;
    }
}
