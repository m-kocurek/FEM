import java.util.ArrayList;
import java.util.List;

public class Grid {
    private List<Node> ND = null;
    private List<Element> EL = null;
    private GlobalData globalData;
    private int nh, ne;

    public Grid(GlobalData globalData) {
        this.globalData = globalData;
        nh = globalData.getNh();
        ne = globalData.getNe();
        ND = new ArrayList<>(nh);
        EL = new ArrayList<>(ne); }

    public void generateGrid(int lpktc){
        double dx = globalData.getDx();
        double dy = globalData.getDy();
        for (int i = 0; i < globalData.getnB(); ++i)
            for (int j = 0; j < globalData.getnH(); ++j){
                double x = i * dx;
                double y = j * dy;
                boolean status = false;
                if(x == 0.00000000 || y == 0.00000 || x == globalData.getB() || y == globalData.getH()){
                    status = true;
                }
                ND.add(new Node(x, y, i * globalData.getnH() + j, status, globalData.getT0() ));
            }

        for (int i = 0 ; i < globalData.getnB() - 1; ++i){
            for (int j = 0; j < globalData.getnH() - 1; ++j){

                int [] tab_ID = new int[4];
                tab_ID[0] = j + i * globalData.getnH();
                tab_ID[3] = tab_ID[0] + 1;
                tab_ID[1] = j + (i+1) * globalData.getnH();
                tab_ID[2] = tab_ID[1] + 1;
                Node [] nodes = new Node[4];
                int z = 0;
                for(int nodeId: tab_ID){
                    nodes[z] = ND.get(nodeId);
                    z++;
                }
                EL.add(new Element(tab_ID, nodes));
            }
        }
    }

    public List getND() {
        return ND;
    }

    public List getEL() {
        return EL;
    }

}
