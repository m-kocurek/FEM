
import java.util.LinkedList;
import java.util.List;
public class Element {

    private int [] IDArray = null;
    private Node [] nodes = null;
    private Sides[] sides = null;
    private int numberOfSidewithBC = 0;
    private List<Integer> IDOfSideswithBC = null;

    public Element(int[] IDArray, Node[] nodes) {
        this.IDArray = IDArray;
        this.nodes = nodes;

        IDOfSideswithBC = new LinkedList<>();
        sides = new Sides[4];
            sides[0] = new Sides(nodes[3], nodes[0]);
            sides[1] = new Sides(nodes[0], nodes[1]);
            sides[2] = new Sides(nodes[1], nodes[2]);
            sides[3] = new Sides(nodes[2], nodes[3]);

            for(Sides side: sides)
                if(side.getSide()[0].isStatus() && side.getSide()[1].isStatus() )
                     numberOfSidewithBC++;

            for(int i = 0; i < 4; i++)
                if(sides[i].getSide()[0].isStatus() && sides[i].getSide()[1].isStatus())
                    IDOfSideswithBC.add(i);
    }

    public Sides getSideOfId(int id){
        if(id >= 0 && id <= sides.length) return sides[id];
        else return null;
    }

    public int getNumberOfSideswithBC() {
        return numberOfSidewithBC;
    }

    public int[] getIDArray() {
        return IDArray;
    }

    public List<Integer> getIDOfSideswithBC() {
        return IDOfSideswithBC;
    }

    public Sides[] getSides() {
        return sides;
    }
    public void setSides(Sides[] side){
        sides = side;
    }


}
