public class Sides {

    private final Node[] side;
    private double [][] N_Fun_Vals;

    public Sides(Node node1, Node node2) {
        this.side = new Node[2];
        this.side[0] = node1;
        this.side[1] = node2;
    }
public Sides(Node node1, Node node2, Node node3){
        this.side=new Node[3];
        this.side[0]=node1;
        this.side[1]=node2;
        this.side[2]=node3;
}

    public Sides(Node node1, Node node2, Node node3, Node node4){
        this.side=new Node[4];
        this.side[0]=node1;
        this.side[1]=node2;
        this.side[2]=node3;
        this.side[3]=node4;
    }

    public void setShapeFunctionVals(double[][] shapeFunctionVals) {
        this.N_Fun_Vals = shapeFunctionVals;
    }

    public double[][] getN_Fun_Vals() {
        return N_Fun_Vals;
    }

    public Node[] getSide() {
        return side;
    }

}
