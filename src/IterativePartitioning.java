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
    ArrayList<double[]> data; //List of input data for cmi computation, currently works for dimension=3
    int pointsSize;  //number of data points
    int dimension;  //dimension of data, must be 3
    int dimCombinat; //2^dimension
    double [] magicNumbers; //Magic number used in the original implementation claimed to be chi2 statistics but they are not
    ArrayList<Integer[]> inverseSortedIndices;  //Instead of double values of points their order in every dimension is used for CMI computation
    ArrayList<Integer[]> sortedIndices;  //indices of points sorted according to their original double value
                        //Usually, only third array from this structure is used.
    int [][] combArray;  //in original code called Imm. This array helps to assign correct bound to new cubes

    LinkedList<Cube> usedCubes;  //list of cubes that were used for CMI computation
    public int maxDeep; //Maximal deep of splitting space into cubes
    TreeNode treeRoot;


    /**
     *
     * @param inputData input data
     * @param toRecompute If toRecompute[i] than it is needed to compute inverseSortedIndices and sortedIndices in the
     *                   i-th coordinate directly from inputData, otherwise it will be reused from other input variables
     * @param oldSorted data structure that can be partially of fully reused for sortedIndices
     * @param oldInverse data structure that can be partially of fully reused for inverseSortedIndices
     */
    public  IterativePartitioning(ArrayList<double[]> inputData,boolean[] toRecompute,ArrayList<Integer[]> oldSorted,ArrayList<Integer[]> oldInverse){
        data=inputData;
        pointsSize=inputData.get(0).length;
        dimension=inputData.size();
        dimCombinat=(int)Math.pow(2,dimension);
        magicNumbers= new double[]{0, 7.81, 13.9, 25.0, 42.0};
        if(toRecompute==null){
            computeInverseOrders();
        }
        else{
            computeInverseOrders2(oldSorted,oldInverse,toRecompute);
        }
        initCombArray();

        usedCubes=new LinkedList<>();
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
        dimension=inverse.size();
        dimCombinat=(int)Math.pow(2,dimension);
        magicNumbers= new double[]{0, 7.81, 13.9, 25.0, 42.0};
        this.sortedIndices=sorted;
        this.inverseSortedIndices=inverse;

        usedCubes=new LinkedList<>();
        maxDeep=Integer.MAX_VALUE;
        initCombArray();
    }


    /**
     * Basic constructor that initializes all necessary structures based on double valued data
     * @param inputData
     */
    public IterativePartitioning(ArrayList<double[]> inputData){
       this(inputData, null, null, null);
    }

    /**
     * Computes sortedIndices and inverseSortedIndices based on input double valued data
     */
    private void computeInverseOrders(){
        inverseSortedIndices = new ArrayList<>();
        sortedIndices=new ArrayList<>();
        for (double[] doubles : data) {
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
    private void computeInverseOrders2(ArrayList<Integer[]> oldSorted,ArrayList<Integer[]> oldInverse,boolean[] toRecompute){
        inverseSortedIndices = new ArrayList<>();
        sortedIndices=new ArrayList<>();
        for (int i = 0; i < dimension; i++) {
            if(toRecompute[i]){
                double [] doubles=data.get(i);
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
    private Cube getInitialCube(){

        int start = 0;
        int end = pointsSize - 1;
        int[] initialCoordinates = {start, end, start, end, start, end};
        Cube initialCube = new Cube(initialCoordinates,0);
        for (int i = 0; i < pointsSize; i++) {
            initialCube.addPoint(i);
        }
        return initialCube;
    }

    private void initCombArray(){
        combArray =new int[dimCombinat][dimension];
        for (int i = 0; i < dimCombinat; i++) {
            int temp = 1;
            for (int j = 0; j < dimension; j++) {
                combArray[i][j] = (i & temp) >> j;
                temp = temp * 2;
            }
        }
    }

    /**
     * Given a cube and coordinates of its center, create eight new cubes.
     * @param currentCube cube to be split
     * @param center center of the cube
     * @return list of eight new cubes
     */
    private ArrayList<Cube> getNewCubes(Cube currentCube,int [] center){

        int[][] currentCoordinates = currentCube.getCoordinates();
        ArrayList<Cube> newCubes=new ArrayList<>();

        for (int i = 0; i < dimCombinat; i++) {
            int [] coordinates=new int[2*dimension];
            for (int j = 0; j < dimension; j++) {
                if(combArray[i][j]>0){
                    coordinates[2*j]=center[j]+1;
                    coordinates[2*j+1]=currentCoordinates[j][1];
                }
                else{
                    coordinates[2*j+1]=center[j];
                    coordinates[2*j]=currentCoordinates[j][0];
                }
            }
            Cube newCube=new Cube(coordinates,currentCube.getDepth()+1);
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
        Cube initialCube=getInitialCube();
        stack.push(initialCube);

        double cmi=0;
        while (!stack.isEmpty()) {
            Cube currentCube=stack.pop();
            ArrayList<Integer> currentPoints = currentCube.getPoints();
            int pointsNumber = currentPoints.size();
            int[] center = currentCube.findCenter();
            ArrayList<Cube> newCubes=getNewCubes(currentCube,center);

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
                        usedCubes.add(smallCube);
                        cmi+=cmiInCube(smallCube);
                    }
                }
            }
            else{
                usedCubes.add(currentCube);
                cmi+=cmiInCube(currentCube);
            }


        }
        cmi=cmi/pointsSize;

        return cmi;

    }


    /**
     * Evaluates whether currentCube or its split: newCubes should be used for CMI computation
     * This method is basically not reliable because it incorporates some strange heuristic based for example on
     * some magic value in magicNumbers[dimension-1] which is claimed to from chi2 statistical table but it is not
     * @param currentCube
     * @param newCubes
     * @return if true, use new Cubes (split), if false, use currentCube (do not split)
     */
    private boolean useSplit(Cube currentCube,ArrayList<Cube> newCubes){  //Can be implemented as abstract and other classes can extend this
        if(currentCube.getDepth()>=maxDeep){
            return false;
        }
        int pointsNumber=currentCube.getPoints().size();     //Or those methods can be a part of an interface and I could use classes implementing the interface
        double expected=(double) pointsNumber /dimCombinat;
        double tst=0;
        for (int i = 0; i < dimCombinat; i++) {
            int size = newCubes.get(i).getPoints().size();
            tst+=(size -expected)*(size -expected);
        }
        tst=dimCombinat*tst/ pointsNumber;
       // System.out.println(tst);
        return (tst>magicNumbers[dimension-1]);
    }

    public boolean pushToStack(Cube newCube){
        return newCube.getPoints().size()>dimCombinat;
    }


    /**
     * for debugging purposes
     * @param point index of point in data and in inverseSortedIndices
     * @return returns coordinates of the point, that is order of given point in every dimension
     */
    public int[] getPointCoordinates(int point){
        int[] coordinates=new int[3];
        coordinates[0]=inverseSortedIndices.get(0)[point];
        coordinates[1]=inverseSortedIndices.get(1)[point];
        coordinates[2]=inverseSortedIndices.get(2)[point];
        return coordinates;
    }


    public double cmiInUniformGrid(int numberOfBins){
        int binSize=pointsSize/numberOfBins;
        int rest=pointsSize%numberOfBins;
        int usedNumberOfBins=numberOfBins;
        if(rest>=0.25*binSize&&rest>=2){
            usedNumberOfBins++;
        }
        TreeMap<Integer,Cube> tmCubes=new TreeMap<>();
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
    public double iterativePartitioningBetter(){
        Stack<Cube> stack = new Stack<>();
        Stack<TreeNode> stack2 = new Stack<>();
        Cube initialCube=getInitialCube();
        TreeNode rootNode=new TreeNode(initialCube);
        stack2.push(rootNode);
        stack.push(initialCube);
        treeRoot=rootNode;

        double cmi=0;
        while (!stack.isEmpty()) {
            Cube currentCube = stack.pop();
            TreeNode currentNode=stack2.pop();
            ArrayList<Integer> currentPoints = currentCube.getPoints();
            int pointsNumber = currentPoints.size();
            int[] center = currentCube.findCenter();
            ArrayList<Cube> newCubes = getNewCubes(currentCube, center);
            for (int i = 0; i < pointsNumber; i++) {
                int sum=0;
                int temp=1;
                Integer pointOrder = currentPoints.get(i);
                for (int j = 0; j < dimension; j++) {
                    if(inverseSortedIndices.get(j)[pointOrder]>center[j]) sum+=temp;
                    temp*=2;
                }
                newCubes.get(sum).addPoint(pointOrder);
            }

            if(useSplit(currentCube,newCubes)||currentCube.equals(initialCube)){
                ArrayList<TreeNode> newNodes=new ArrayList<>();
                for (int i = 0; i < dimCombinat; i++) {
                    TreeNode tn=new TreeNode(newCubes.get(i));
                    newNodes.add(tn);
                }
                //TreeNode currentNode=stack2.pop();
                currentNode.reinit(center,newNodes);
                for (int i = 0; i < dimCombinat; i++) {

                    Cube smallCube = newCubes.get(i);
                    if (pushToStack(smallCube)){
                        stack.push(smallCube);
                        stack2.push(newNodes.get(i));
                    }
                    else if(!smallCube.getPoints().isEmpty()){
                        usedCubes.add(smallCube);
                        //cmi+=cmiInCube(smallCube);
                        cmi+=cmiInCubeBetter(smallCube);
                    }
                    else if(smallCube.getPoints().isEmpty()){
                        usedCubes.add(smallCube);
                    }
                }
            }
            else{
                //cmi+=cmiInCube(currentCube);
                cmi+=cmiInCubeBetter(currentCube);
                usedCubes.add(currentCube);
            }


        }
        cmi=cmi/pointsSize;

        return cmi;
    }


    /**
     * Computes CMI in one cube (grid cell) in the same way as in Cupid. The aim was to implement it such that the
     * time complexity is close to original implementation as much as possible
     * @param cube cube in which we compute CMI
     * @return a value that needs to be divided by pointsSize to obtain the contribution to final CMI value
     */
    public double cmiInCube(Cube cube){

        int[][] coordinates = cube.getCoordinates();
        int sizePointsInCube = cube.getPoints().size();

        int nz=coordinates[2][1]-coordinates[2][0]+1;
        int nxz=0;
        int nyz=0;
        double cmi=0;
        for (int i = 0; i < pointsSize; i++) {
            int x=inverseSortedIndices.get(0)[i];
            int y=inverseSortedIndices.get(1)[i];
            int z=inverseSortedIndices.get(2)[i];
            if(x>=coordinates[0][0]&&x<=coordinates[0][1]){
                if(z>=coordinates[2][0]&&z<=coordinates[2][1]){
                    nxz++;
                }
            }
            if(y>=coordinates[1][0]&&y<=coordinates[1][1]){
                if(z>=coordinates[2][0]&&z<=coordinates[2][1]){
                    nyz++;
                }
            }
        }

        double cond=1;
        if(nxz*nyz>0&&nz*sizePointsInCube>0){

            cond=((double)nz* sizePointsInCube)/(nxz*nyz);
        }
        cmi=sizePointsInCube*Math.log(cond);

        return cmi;
    }


    /**
     * This is a better version of CMI computation in one cube (grid cell)
     * @param cube
     * @return a value that needs to be divided by pointsSize to obtain the contribution to final CMI value
     */
    public double cmiInCubeBetter(Cube cube){

        int[][] coordinates = cube.getCoordinates();
        int sizePointsInCube = cube.getPoints().size();
        int nz=coordinates[2][1]-coordinates[2][0]+1;
        double cond=1;
        double cmi = 0;

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

         if(nxz*nyz>0&&nz*sizePointsInCube>0){

                cond = ((double) nz * sizePointsInCube) / (nxz * nyz);  //In a better version the condition that tests against zero could be used earlier to prevent useless computation
            }

        cmi=sizePointsInCube*Math.log(cond);

        return cmi;
    }

    /**
     * For debugging purposes. Tests whether all points assigned to a cube are within its bounds
     * @param cube
     */
    public void testCube(Cube cube){
        boolean agree=true;
        int[][] coordinates = cube.getCoordinates();
        for (Integer point : cube.getPoints()) {
            for (int i = 0; i < dimension; i++) {
                int currentCoordinate=inverseSortedIndices.get(i)[point];
                if(currentCoordinate<coordinates[i][0]||currentCoordinate>coordinates[i][1]){
                    System.out.println("error in point split");
                    System.out.println("Point "+point+" "+i+"-th coordinate is "+currentCoordinate);
                    System.out.println("but it has to be between "+coordinates[i][0]+", "+coordinates[i][1]);
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
//        HashSet<Integer> usedPoints=new HashSet<>();
//        for (Cube cube : cubes) {
//            testCube(cube);
//            usedPoints.addAll(cube.getPoints());
//        }
//        if(usedPoints.size()!=pointsSize){
//            System.out.println("Not all points assigned");
//        }
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
//        Stack<TreeNode> stack=new Stack<>();
//        stack.push(treeRoot);
//        while(!stack.isEmpty()){
//            TreeNode tn=stack.pop();
//            boolean isLeaf=tn.isLeaf;
//            if(!isLeaf){
//                for (TreeNode node : tn.nextNodes) {
//                    stack.push(node);
//                }
//            }
//            else{
//                tn.cube.getPoints().clear();
//            }
//        }
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
//        stack=new Stack<>();
//        stack.push(treeRoot);
//        while(!stack.isEmpty()){
//            TreeNode tn=stack.pop();
//            boolean isLeaf=tn.isLeaf;
//            if(!isLeaf){
//                for (TreeNode node : tn.nextNodes) {
//                    stack.push(node);
//                }
//            }
//            else{
//                //testCube(tn.cube);
//                cmi+=cmiInCubeBetter(tn.cube);
//            }
//        }
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



}
