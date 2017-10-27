package org.mpii.jami.cmi;

import java.util.ArrayList;

/**
 * Created by fuksova on 1/28/16.
 * Structure used in for computation of CMI in a fixed grid. Should be replaced by another approach because it is not
 * effective.
 */
public class TreeNode {
    int[] cubeCenter;
    ArrayList<TreeNode> nextNodes;
    boolean isLeaf;
    Cube cube;



    public TreeNode(Cube cube){
        this.cube=cube;
        cubeCenter=null;
        isLeaf=true;
        cubeCenter=null;
    }

    public void reinit(int[] cubeCenter,ArrayList<TreeNode> cubes){
        isLeaf=false;
        cube=null;
        this.cubeCenter=cubeCenter;
        nextNodes =cubes;
    }
}
