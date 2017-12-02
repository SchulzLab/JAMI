package org.mpii.jami.cmi;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.mpii.jami.helpers.AdditionalComputations;
import org.mpii.jami.helpers.ComparatorForIndices;

import java.util.*;

/**
 * Created by fuksova on 11/9/15.
 * Iterative partitioning computed in the same way as in Cupid. One of its drawbacks is that it assigns nonzero
 * CMI value when all values in data.get(1) are equal to zero. In combination with other principles used in Cupid,
 * such data sample can be even evaluated as significant.
 * For this computation original double values of points are replaced by their integer order in every dimension. This
 * information is stored in inverseSortedIndices
 * Points are stored in CMI grid cells in form of their indices in the original input data structure.
 */
public class IterativePartitioning {
    private static final Logger logger = LogManager.getLogger("JAMI");

    private ArrayList<List<Double>> data; //List of input data for cmi computation, currently works for dimension=3
    private int pointsSize;  //number of data points
    private double [] chiSquareCutoffs;
    private boolean considerZeros;

    public ArrayList<Integer[]> getInverseSortedIndices() {
        return inverseSortedIndices;
    }

    public ArrayList<Integer[]> getSortedIndices() {
        return sortedIndices;
    }

    private ArrayList<Integer[]> inverseSortedIndices;  //Instead of double values of points their order in every dimension is used for CMI computation
    private ArrayList<Integer[]> sortedIndices;  //indices of points sorted according to their original double value
                        //Usually, only third array from this structure is used.
    private HashMap<Integer, int [][]> combArrays = new HashMap<Integer, int[][]>();  //in original code called Imm. This array helps to assign correct bound to new cubes

    private int maxDeep; //Maximal deep of splitting space into cubes
    private ArrayList<HashSet<Integer>> minimumValues;


    public boolean ConsiderZeros() {
        return considerZeros;
    }

    public void setConsiderZeros(boolean considerZeros) {
        this.considerZeros = considerZeros;
    }

    public int getMaxDeep() {
        return maxDeep;
    }

    public void setMaxDeep(int maxDeep) {
        this.maxDeep = maxDeep;
    }
    /**
     *
     * @param inputData input data
     * @param toRecompute If toRecompute[i] than it is needed to compute inverseSortedIndices and sortedIndices in the
     *                   i-th coordinate directly from inputData, otherwise it will be reused from other input variables
     * @param oldSorted data structure that can be partially of fully reused for sortedIndices
     * @param oldInverse data structure that can be partially of fully reused for inverseSortedIndices
     */
    public  IterativePartitioning(ArrayList<List<Double>> inputData,boolean[] toRecompute,ArrayList<Integer[]> oldSorted,ArrayList<Integer[]> oldInverse){
        data=inputData;
        pointsSize=inputData.get(0).size();
        chiSquareCutoffs = new double[]{0, 7.81, 13.9, 25.0, 42.0};
        if(toRecompute==null){
            computeInverseOrders();
        }
        else{
            computeInverseOrders2(oldSorted,oldInverse,toRecompute, inputData.size());
        }
        this.combArrays.put(2, initCombArray(2, 4));
        this.combArrays.put(3, initCombArray(3, 8));

        maxDeep=Integer.MAX_VALUE;
    }

    /**
     * No double data are need. The whole computation is performed on inverseSortedIndices
     * @param sorted initial value for sortedIndices
     * @param inverse initial value for inverseSortedIndices
     */
    public IterativePartitioning(ArrayList<Integer[]> sorted,ArrayList<Integer[]> inverse){
        data=null;
        pointsSize=inverse.get(0).length;
        chiSquareCutoffs = new double[]{0, 7.81, 13.9, 25.0, 42.0};
        this.sortedIndices=sorted;
        this.inverseSortedIndices=inverse;

        maxDeep=Integer.MAX_VALUE;

        this.combArrays.put(2, initCombArray(2, 4));
        this.combArrays.put(3, initCombArray(3, 8));
    }

    /**
     * Basic constructor that initializes all necessary structures based on double valued data
     * @param inputData
     */
    public IterativePartitioning(ArrayList<List<Double>> inputData){
       this(inputData, null, null, null);
    }

    /**
     * Computes sortedIndices and inverseSortedIndices based on input double valued data
     */
    private void computeInverseOrders(){
        inverseSortedIndices = new ArrayList<>();
        sortedIndices=new ArrayList<>();
        for (List<Double> doubles : data) {
            ComparatorForIndices comp = new ComparatorForIndices(doubles);
            Integer[] arrayToSort = comp.getArrayToSort();
            Arrays.sort(arrayToSort, comp);
            sortedIndices.add(arrayToSort);
            Integer[] inverse = comp.getInverseArray(arrayToSort);
            inverseSortedIndices.add(inverse);
        }
    }


    /**
     * Used when part of inverseSortedIndices and sortedIndices can be reused from previously computed values
     *
     * @param oldSorted partially of fully used for sortedIndices computation
     * @param oldInverse partially of fully used for inverseSortedIndices computation
     * @param toRecompute toRecompute[i]==true means that i-th coordinate data for inverseSortedIndices and sortedIndices
     *                    should be recomputed. Otherwise, values from oldInverse and oldSorted should be used.
     */
    private void computeInverseOrders2(ArrayList<Integer[]> oldSorted,
                                       ArrayList<Integer[]> oldInverse,boolean[] toRecompute, int dimension){
        inverseSortedIndices = new ArrayList<>();
        sortedIndices=new ArrayList<>();
        for (int i = 0; i < dimension; i++) {
            if(toRecompute[i]){
                List<Double> doubles=data.get(i);
                ComparatorForIndices comp = new ComparatorForIndices(doubles);
                Integer[] arrayToSort = comp.getArrayToSort();
                Arrays.sort(arrayToSort, comp);
                sortedIndices.add(arrayToSort);
                Integer[] inverse = comp.getInverseArray(arrayToSort);
                inverseSortedIndices.add(inverse);
            }
            else{
                sortedIndices.add(oldSorted.get(i));
                inverseSortedIndices.add(oldInverse.get(i));
            }
        }

    }

    /**
     *
     * @return Initial cube that has bounds 1..pointsSize in every dimension and contains all input points.
     */
    public static Cube getInitialCube(int pointsSize){

        int start = 0;
        int end = pointsSize - 1;
        int[] initialCoordinates = {start, end, start, end, start, end};
        Cube initialCube = new Cube(initialCoordinates,0);
        for (int i = 0; i < pointsSize; i++) {
            initialCube.addPoint(i);
        }
        return initialCube;
    }

    private int[][] initCombArray(int dimension, int dimCombinat) {
        int[][] combArray = new int[dimCombinat][dimension];
        for (int i = 0; i < dimCombinat; i++) {
            int temp = 1;
            for (int j = 0; j < dimension; j++) {
                combArray[i][j] = (i & temp) >> j;
                temp = temp * 2;
            }
        }
        return combArray;
    }

    /**
     * Given a cube and coordinates of its center, create 2^dim new hypercubes.
     * @param currentCube cube to be split
     * @param center center of the cube
     * @return list of 2^dim new cubes
     */
    private ArrayList<Cube> getNewCubes(Cube currentCube,int [] center){

        int[][] currentCoordinates = currentCube.getCoordinates();
        ArrayList<Cube> newCubes=new ArrayList<>();
        int dimension = currentCube.getDimensions().length;

        for (int i = 0; i < (int) Math.pow(2, dimension); i++) {
            int [] coordinates=new int[2*dimension];
            for (int j = 0; j < dimension; j++) {
                if(combArrays.get(dimension)[i][j]>0){
                    coordinates[2*j]=center[j]+1;
                    coordinates[2*j+1]=currentCoordinates[j][1];
                }
                else{
                    coordinates[2*j+1]=center[j];
                    coordinates[2*j]=currentCoordinates[j][0];
                }
            }
            Cube newCube=new Cube(coordinates,currentCube.getDepth()+1);
            newCube.setDimensions(currentCube.getDimensions());
            newCube.setPointsAfterZeroSplit(currentCube.getPointsAfterZeroSplit());
            newCubes.add(newCube);
        }
        return newCubes;
    }


    /**
     * This implementation of iterative partitioning tries to copy as accurate as possible the original implementation
     * in Cupid in the sense of time complexity.
     * @return CMI value
     */
    public double naivePartitioning() {
        Stack<Cube> stack = new Stack<>();
        Cube initialCube=getInitialCube(this.pointsSize);
        stack.push(initialCube);
        int dimCombinat = (int) Math.pow(2, initialCube.getDimensions().length);

        double cmi=0;
        while (!stack.isEmpty()) {
            Cube currentCube=stack.pop();
            ArrayList<Integer> currentPoints = currentCube.getPoints();
            int pointsNumber = currentPoints.size();
            int[] center = currentCube.findCenter();
            ArrayList<Cube> newCubes=getNewCubes(currentCube,center);
            int dimension = currentCube.getDimensions().length;
            int[][] J = new int[pointsNumber][dimension];
            for (int i = 0; i < dimension; i++) {
                Integer[] currentColumn = inverseSortedIndices.get(i);
                for (int j = 0; j < pointsNumber; j++) {
                    if(currentColumn[currentPoints.get(j)]>center[i]) J[j][i]=1;
                }
            }
            for (int d = 0; d < dimCombinat; d++) {
                for (int i = 0; i < pointsNumber; i++) {
                    int index=0;
                    int temp=1;
                    for (int j = 0; j < dimension; j++) {
                        index += J[i][j] * temp;
                        temp *= 2;
                    }
                    if(index==d){   //this is and artificially added condition to ensure similar complexity. Normally,
                                      //it could be done without the condition and without the outer cycle
                        newCubes.get(index).addPoint(currentPoints.get(i));
                    }
                }

            }

            for (Cube newCube : newCubes) {
                testCube(newCube);
            }

            if(useSplit(currentCube,newCubes)||currentCube.equals(initialCube)){
                for (int i = 0; i < dimCombinat; i++) {
                    Cube smallCube = newCubes.get(i);
                    if (pushToStack(smallCube)){
                        stack.push(smallCube);
                    }
                    else if(!smallCube.getPoints().isEmpty()){
                        cmi+=cmiInCube(smallCube);
                    }
                }
            }
            else{
                cmi+=cmiInCube(currentCube);
            }


        }
        cmi=cmi/pointsSize;

        return cmi;

    }


    /**
     * Evaluates whether currentCube or its split: newCubes should be used for CMI computation
     * This method is basically not reliable because it incorporates some strange heuristic based for example on
     * some magic value in chiSquareCutoffs[dimension-1] which is claimed to from chi2 statistical table but it is not
     * @param currentCube
     * @param newCubes
     * @return if true, use new Cubes (split), if false, use currentCube (do not split)
     */
    private boolean useSplit(Cube currentCube,ArrayList<Cube> newCubes){
        if(currentCube.getDepth()>=maxDeep){
            return false;
        }
        int dimCombinat = (int) Math.pow(2, currentCube.getDimensions().length);
        int pointsNumber=currentCube.getPoints().size();
        double expected=(double) pointsNumber / dimCombinat;
        double tst=0;
        for (int i = 0; i < dimCombinat; i++) {
            int size = newCubes.get(i).getPoints().size();
            tst+=(size -expected)*(size -expected);
        }
        tst=dimCombinat*tst/ pointsNumber;

        return (tst> chiSquareCutoffs[currentCube.getDimension()-1]);
    }

    public boolean pushToStack(Cube newCube){
        return newCube.getPoints().size()>Math.pow(2, newCube.getDimension());
    }

    public double cmiInUniformGrid(int numberOfBins){
        int binSize=pointsSize/numberOfBins;
        int rest=pointsSize%numberOfBins;
        int usedNumberOfBins=numberOfBins;
        int dimension = data.size();
        if(rest>=0.25*binSize&&rest>=2){
            usedNumberOfBins++;
        }
        int[] coeffs=new int[dimension];

        coeffs[0]=usedNumberOfBins*usedNumberOfBins;
        coeffs[1]=usedNumberOfBins;
        coeffs[2]=1;
        int numCubes=usedNumberOfBins*usedNumberOfBins*usedNumberOfBins;
        Cube [] listCubes=new Cube[numCubes];
        int counter=0;
        int[] keys1=new int[dimension];
        for (int i = 0; i <usedNumberOfBins; i++) {
            keys1[0]=i;
            for (int j = 0; j < usedNumberOfBins; j++) {
                keys1[1]=j;
                for (int k = 0; k < usedNumberOfBins; k++) {
                    keys1[2]=k;
                    int[] coordinates=new int[2*dimension];
                    for (int l = 0; l < dimension; l++) {
                        coordinates[2*l]=keys1[l]*binSize;
                        coordinates[2*l+1]=(keys1[l]+1)*binSize-1;
                    }
                    Cube cube=new Cube(coordinates,0);
                    listCubes[counter]=cube;
                    counter++;

                }
            }
        }
        for (int i = 0; i < pointsSize; i++) {
            int[] keys=new int[dimension];
            Integer hash=0;
            for (int j = 0; j < dimension; j++) {
                keys[j]=Math.min(inverseSortedIndices.get(j)[i]/binSize,usedNumberOfBins-1);
                hash+=keys[j]*coeffs[j];
            }
            Cube cube=listCubes[hash];

            cube.addPoint(i);

        }
        double cmi=0;

        for (int i = 0; i < listCubes.length; i++) {
            cmi+=cmiInCubeBetter(listCubes[i]);
        }
        return cmi/pointsSize;

    }


    /**
     * The same principle of iterative partitioning as in Cupid but with some implementation enhancements
     * @return CMI value
     */
    public double iterativePartitioningBetter(Cube initialCube){

        Stack<Cube> stack = new Stack<>();

        if(this.considerZeros) stack = splitCubeWithZeros(initialCube);
        else stack.push(initialCube);

        double cmi=0;
        while (!stack.isEmpty()) {
            Cube currentCube = stack.pop();
            ArrayList<Integer> currentPoints = currentCube.getPoints();
            int pointsNumber = currentPoints.size();
            int[] center = currentCube.findCenter();
            ArrayList<Cube> newCubes = getNewCubes(currentCube, center);
            for (int i = 0; i < pointsNumber; i++) {
                int sum = 0;
                int temp = 1;
                Integer pointOrder = currentPoints.get(i);
                for(int j = 0; j < currentCube.getDimensions().length; j++) {
                    if(inverseSortedIndices.get(currentCube.getDimensions()[j])[pointOrder]>center[j]) sum+=temp;
                    temp*=2;
                }
                newCubes.get(sum).addPoint(pointOrder);
            }

            if(useSplit(currentCube,newCubes)||currentCube.equals(initialCube)){
                for (int i = 0; i < (int) Math.pow(2, currentCube.getDimensions().length); i++) {

                    Cube smallCube = newCubes.get(i);
                    if (pushToStack(smallCube)){
                        stack.push(smallCube);
                    }
                    else if(!smallCube.getPoints().isEmpty()){
                        cmi+=cmiInCubeBetter(smallCube);
                    }
                }
            }
            else{
                cmi+=cmiInCubeBetter(currentCube);
            }


        }
        cmi=cmi/pointsSize;

        return cmi;
    }

    private Stack<Cube> splitCubeWithZeros(Cube initialCube) {

        Stack<Cube> cubesToConsiderForSplitting = new Stack<>();
        cubesToConsiderForSplitting.push(initialCube);
        Stack<Cube> finalCubes = new Stack<>();
        ArrayList<HashSet<Integer>> minimumValues = new ArrayList<>();

        if(this.data != null) {
            for (int dim = 0; dim < initialCube.getDimension(); dim++) {
                minimumValues.add(AdditionalComputations.ranksOfMinimumValues(this.data.get(dim)));
            }
            setMinimumValues(minimumValues);
        } else{
            minimumValues = getMinimumValues();
        }

        while(!cubesToConsiderForSplitting.empty()){

            Cube currentCube = cubesToConsiderForSplitting.pop();

            if(currentCube.getDimension() == 1) continue; //currently not handling one dimensional cubes

            boolean noZeros = true;

            //in each dimension we split off a hypercube with zeros which we will flatten to dim - 1 and process further
            for(int dim : currentCube.getDimensions()) {

                if(dim == currentCube.getDimension()) break; // dimension no longer relevant

                //check which of the points in the current cube are at the minimum value
                HashSet<Integer> zeroPointsInCube = new HashSet<>(minimumValues.get(dim));
                zeroPointsInCube.retainAll(currentCube.getPoints());

                if (zeroPointsInCube.size() > 1) { //one value at minimum means no duplicates
                    if(currentCube.getDimension() == 1) break; //only zeros in this cube, cmi is 0.

                    noZeros = false;
                    ArrayList<Integer> nonZeroPoints = new ArrayList<>(currentCube.getPoints());
                    nonZeroPoints.removeAll(zeroPointsInCube);

                    int[][] currentCoordinates = currentCube.cloneCoordinates();
                    int[][] coordinatesOfZeroCube = new int[currentCube.getDimension() -1][2];

                    // we need to figure out which dimension to ignore in building the square
                    int newDim = 0;
                    int[] dimensions = new int[currentCube.getDimension() - 1];
                    for (int oldDim = 0; oldDim < currentCube.getDimension(); oldDim++) {
                        if (oldDim != dim) {
                            coordinatesOfZeroCube[newDim][0] = currentCoordinates[oldDim][0];
                            coordinatesOfZeroCube[newDim][1] = currentCoordinates[oldDim][1];
                            dimensions[newDim++] = oldDim;
                        }
                    }

                    currentCoordinates[dim][0] = zeroPointsInCube.size();

                    Cube nonZeroCube = new Cube(currentCoordinates, currentCube.getDepth() + 1, currentCube.getDimensions());
                    Cube zeroCube = new Cube(coordinatesOfZeroCube, currentCube.getDepth() + 1, dimensions);

                    zeroCube.addPoints(zeroPointsInCube);
                    zeroCube.setPointsAfterZeroSplit(zeroPointsInCube);

                    nonZeroCube.addPoints(nonZeroPoints);
                    if(nonZeroCube.getDimension() < 3){
                        HashSet<Integer> lastZeros = new HashSet<>(currentCube.getPointsAfterZeroSplit());
                        lastZeros.retainAll(nonZeroCube.getPoints());
                        nonZeroCube.setPointsAfterZeroSplit(lastZeros);
                    }
                    cubesToConsiderForSplitting.push(zeroCube);

                    if (nonZeroPoints.size() > 0) currentCube = nonZeroCube;
                    else break;
                }
            }

            if(noZeros) finalCubes.push(currentCube);
            else if(currentCube.getPoints().size() >0) cubesToConsiderForSplitting.push(currentCube);
        }

        return finalCubes;
    }


    /**
     * Computes CMI in one cube (grid cell) in the same way as in Cupid. The aim was to implement it such that the
     * time complexity is close to original implementation as much as possible
     *
     * @param cube cube in which we compute CMI
     * @return a value that needs to be divided by pointsSize to obtain the contribution to final CMI value
     */
    public double cmiInCube(Cube cube) {

        if (cube.getDimensions().length < 3) return (cmiInSquare(cube));

        int[][] coordinates = cube.getCoordinates();
        int sizePointsInCube = cube.getPoints().size();

        int nz = coordinates[2][1] - coordinates[2][0] + 1;
        int nxz = 0;
        int nyz = 0;
        double cmi;
        for (int i = 0; i < pointsSize; i++) {
            int x = inverseSortedIndices.get(0)[i];
            int y = inverseSortedIndices.get(1)[i];
            int z = inverseSortedIndices.get(2)[i];
            if (x >= coordinates[0][0] && x <= coordinates[0][1]) {
                if (z >= coordinates[2][0] && z <= coordinates[2][1]) {
                    nxz++;
                }
            }
            if (y >= coordinates[1][0] && y <= coordinates[1][1]) {
                if (z >= coordinates[2][0] && z <= coordinates[2][1]) {
                    nyz++;
                }
            }
        }

        double cond = 1;
        if (nxz * nyz > 0 && nz * sizePointsInCube > 0) {

            cond = ((double) nz * sizePointsInCube) / (nxz * nyz);
        }
        cmi = sizePointsInCube * Math.log(cond);

        return cmi;
    }

    /**
     * This function computes the mutual information on a cube, ignoring one of the dimensions.
     * @param square
     * @return double mutual information of the cube with one dimension equal to zero (no entropy)
     */
    public double cmiInSquare(Cube square){

        int sizePointsInCube = square.getPoints().size();
        if(sizePointsInCube == 0) return(0);

        int[][] coordinates = square.getCoordinates();

        //find out which of the three dimensions to ignore
        int firstDim = square.getDimensions()[0];
        int secondDim = square.getDimensions()[1];
        int flattenedDim = 3 - firstDim - secondDim;

        int nxz = 0;
        int nyz = 0;

        double cond = 1;

        //case 1 and 2, x or y have been flattened. the way we compute the new coordinates we can be sure that the
        //second dimension is always Z in I(X;Y|Z)
        if(flattenedDim < 2){
            for (int i = coordinates[1][0]; i <= coordinates[1][1]; i++) {
                int point = (sortedIndices.get(secondDim)[i]);
                int x = inverseSortedIndices.get(firstDim)[point];
                int y = inverseSortedIndices.get(flattenedDim)[point];
                if(square.getPointsAfterZeroSplit().contains(y)){
                    nyz++;
                }
                if (x >= coordinates[0][0] && x <= coordinates[0][1]) {
                    nxz++;
                }
            }
            if(nxz > 0)
                cond = ((double) sizePointsInCube) * (coordinates[1][1] - coordinates[1][0] + 1) /
                        (nxz * nyz);
        }

        //case 3, z has been flattened. we consider all points that were included in the initial square since
        // all those points lie on the z-axis even after further splits.
        else{
            for (int i : square.getPointsAfterZeroSplit()) {
                int x=inverseSortedIndices.get(firstDim)[i];
                int y=inverseSortedIndices.get(secondDim)[i];

                if(x>=coordinates[0][0]&&x<=coordinates[0][1]){
                    nxz++;
                }
                if(y>=coordinates[1][0]&&y<=coordinates[1][1]){
                    nyz++;
                }
            }
            if(nxz*nyz > 0)
                cond = ((double) sizePointsInCube) * ((double) square.getPointsAfterZeroSplit().size()) /(nxz*nyz);
        }
        double mi = sizePointsInCube*Math.log(cond);

        return mi;
    }

    /**
     * This is a better version of CMI computation in one cube (grid cell)
     * @param cube
     * @return a value that needs to be divided by pointsSize to obtain the contribution to final CMI value
     */
    public double cmiInCubeBetter(Cube cube){

        if (cube.getDimension() < 3) return(cmiInSquare(cube));

        int[][] coordinates = cube.getCoordinates();
        int sizePointsInCube = cube.getPoints().size();
        int nz=coordinates[2][1]-coordinates[2][0]+1;
        if(nz*sizePointsInCube == 0) return(0);

        double cond=1;
        double cmi;

            int nxz = 0;
            int nyz = 0;
            Integer[] xCoords = inverseSortedIndices.get(0);
            Integer[] yCoords = inverseSortedIndices.get(1);
            Integer[] zSorted = sortedIndices.get(2);
            for (int i = coordinates[2][0]; i <= coordinates[2][1]; i++) {
                int point = (zSorted[i]);
                int x = xCoords[point];
                int y = yCoords[point];

                if (x >= coordinates[0][0] && x <= coordinates[0][1]) {
                    nxz++;
                }
                if (y >= coordinates[1][0] && y <= coordinates[1][1]) {
                    nyz++;
                }
            }

         if(nxz*nyz > 0){

                cond = ((double) nz * sizePointsInCube) / (nxz * nyz);
            }

        cmi=sizePointsInCube*Math.log(cond);

        return cmi;
    }

    /**
     * For debugging purposes. Tests whether all points assigned to a cube are within its bounds
     * @param cube
     */
    public void testCube(Cube cube){
        int[][] coordinates = cube.getCoordinates();
        for (Integer point : cube.getPoints()) {
            for (int i = 0; i < cube.getDimensions().length; i++) {
                int currentCoordinate=inverseSortedIndices.get(i)[point];
                if(currentCoordinate<coordinates[i][0]||currentCoordinate>coordinates[i][1]){
                    logger.debug("Error in point split: " +
                    "Point "+point+" "+i+"-th coordinate is "+currentCoordinate +
                            "but it has to be between "+coordinates[i][0]+", "+coordinates[i][1]);
                }
            }
        }
    }

    /**
     * This CMI computation is based on principle that a fixed grid obtained on randomized data can be reused for
     * all other randomization of the same data. The old grid is given in the form of list of cubes. Those cubes can be
     * obtained during normal iterative partitioning procedure which stores all cubes from the final grid into usedCubes.
     * This method needs to be completely reimplemented to be more effective!
     * @param cubes
     * @return
     */
    public double computeCMIInGrid(LinkedList<Cube> cubes){
        if(cubes.isEmpty()){
            throw new IllegalArgumentException("Grid for CMI computation can not be empty");
        }
        double cmi=0;
        for (Cube cube : cubes) {
            cube.getPoints().clear();
            HashSet<Integer> newPointsInCube = selectPointsInCube(cube);
            if(!newPointsInCube.isEmpty()){
                cube.addPoints(newPointsInCube);
                cmi+=cmiInCube(cube);
            }

            //maybe bad idea, x and y or such coordinates can be reused in case of randomized one coordinate
            //clever association of points into correct cube then compute cmi in cube
        }
        cmi=cmi/pointsSize;
        return cmi;
    }

    /**
     * Another approach to computation of CMI in a fixed grid. Needs to be reimplemented.
     * @param treeRoot
     * @param cubes
     * @return
     */
    public double computeCMIInGrid(TreeNode treeRoot,LinkedList<Cube> cubes){
        double cmi=0;
        int dimension = data.size();

        for (Cube cube : cubes) {
            cube.getPoints().clear();
        }
        for (int p = 0; p < pointsSize; p++) {
            int[] coordinates=new int[dimension];
            for (int i = 0; i < dimension; i++) {
                coordinates[i]=inverseSortedIndices.get(i)[p];
            }
            TreeNode currentNode=treeRoot;
            while(!currentNode.isLeaf){
                int[] center = currentNode.cubeCenter;
                int sum=0;
                int temp=1;
                for (int j = 0; j < dimension; j++) {
                    if(coordinates[j]>center[j]) sum+=temp;
                    temp*=2;
                }
                TreeNode nextNode = currentNode.nextNodes.get(sum);
                currentNode=nextNode;
            }
            currentNode.cube.addPoint(p);
            testCube(currentNode.cube);
        }

        for (Cube cube : cubes) {
            cmi+=cmiInCubeBetter(cube);

        }
        cmi=cmi/pointsSize;
        return cmi;
    }



    /**
     * given bounds of cube, select points that belong to it. Used in computeCMIInGrid
     * @param cube Cube to which points should be assigned. List of points in this cube should be empty before calling
     *             this function
     * @return set of indices of points that should be assigned to the cube.
     */
    private HashSet<Integer> selectPointsInCube(Cube cube){  //Use only for distributing points into predefined grid
        int[][] coordinates = cube.getCoordinates();
        HashSet<Integer> set1=new HashSet<>();


        Integer[] xCoords = inverseSortedIndices.get(0);
        Integer[] yCoords = inverseSortedIndices.get(1);
        Integer[] zSorted = sortedIndices.get(2);
        for (int i = coordinates[2][0]; i <= coordinates[2][1]; i++) {
            int point = (zSorted[i]);
            int x = xCoords[point];
            int y = yCoords[point];

            if ((x >= coordinates[0][0] && x <= coordinates[0][1])&&(y >= coordinates[1][0])&&(y<=coordinates[1][1])) {
               set1.add(point);
            }
        }


        return set1;
    }

    public void setMinimumValues(ArrayList<HashSet<Integer>> minimumValues) {
        this.minimumValues = minimumValues;
    }

    public ArrayList<HashSet<Integer>> getMinimumValues() {
        return minimumValues;
    }
}
