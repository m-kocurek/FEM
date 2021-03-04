import com.google.gson.Gson;
import java.lang.*;
import java.io.FileReader;
import java.util.Arrays;
import Jama.Matrix;

public class GlobalData {

    private double H;
    private double B;
    private double dy;
    private double dx;
    private double k;
    private int lpktc;
    private double c;
    private double ro;
    private double alfa;
    private double tau;
    private double t0;
    private double t_alfa;
    private double dtau;
    private int nH;
    private int nB;
    private int nh;
    private int ne;
    private Matrix hLocal;
    private Matrix hGlobal;
    private Matrix cLocal;
    private Matrix cGlobal;
    private Matrix trans_matrix;
    private Matrix Hbc;
    private Matrix shapeFunctionsEta;
    private Matrix shapeFunctionsKsi;
    private Matrix shapeFunctions;
    private Element localElement;
    private Matrix border_cond;
    private Matrix p_local;
    private Matrix p_local_suma;
    private Matrix pGlobal;
    private double det_J_global;
    private Grid grid;
    private double Temp =0.0;

    public GlobalData(  boolean x ){
        GlobalData Data=this.generate_data() ;
        if(Data != null) {
            this.B = Data.getB();
            this.H = Data.getH();
            this.nH = Data.getnH();
            this.nB = Data.getnB();
            this.lpktc=Data.getLpktc();
            this.alfa=Data.getAlfa();
            this.setNh(this.getnH() * this.getnB());
            this.setNe((this.getnB() - 1) * (this.getnH() - 1));
            this.setDx(this.getB() / (this.getnB() - 1));
            this.setDy(this.getH() / (this.getnH() - 1));
            this.hLocal = null;
            this.hGlobal = new Matrix(this.nh, this.nh);
            this.cLocal = null;
            this.cGlobal = new Matrix(this.nh, this.nh);
            this.trans_matrix= new Matrix(4,lpktc*lpktc);
            this.k=Data.getK();
            this.c= Data.getC();
            this.ro = Data.getRo();
            this.alfa= Data.getAlfa();
            this.t0 = Data.getT0();
             this.t_alfa=Data.getT_alfa();
             this.dtau=Data.getDt();
             this.tau=Data.getTau();
            generateShapeFunMatrix(lpktc);
            generateLocalElement();
            grid = new Grid(this);
            grid.generateGrid(lpktc);
        }
    }


    public void generateShapeFunMatrix(int lpktc){

        Point[] points = IntegralPoints.getIntegralPoints( lpktc);

        shapeFunctionsEta = new Matrix(lpktc*lpktc, 4);
        for(int i = 0; i < lpktc*lpktc; i++) {
            shapeFunctionsEta.set(i, 0,  shapeFunction1Eta(points[i].getX()));
            shapeFunctionsEta.set(i, 1,  shapeFunction2Eta(points[i].getX()));
            shapeFunctionsEta.set(i, 2,  shapeFunction3Eta(points[i].getX()));
            shapeFunctionsEta.set(i, 3,  shapeFunction4Eta(points[i].getX()));
        }

        shapeFunctionsKsi = new Matrix(lpktc*lpktc, 4);
        for(int i = 0; i < lpktc*lpktc; i++){
            shapeFunctionsKsi.set(i,0,  shapeFunction1Ksi(points[i].getY()));
            shapeFunctionsKsi.set(i,1,  shapeFunction2Ksi(points[i].getY()));
            shapeFunctionsKsi.set(i,2,  shapeFunction3Ksi(points[i].getY()));
            shapeFunctionsKsi.set(i,3,  shapeFunction4Ksi(points[i].getY()));
        }

        shapeFunctions = new Matrix(lpktc*lpktc,4);
        for(int i = 0; i < lpktc*lpktc; i++){
            shapeFunctions.set(i,0,  shapeFunction1(points[i].getX(), points[i].getY()));
            shapeFunctions.set(i,1,  shapeFunction2(points[i].getX(), points[i].getY()));
            shapeFunctions.set(i,2,  shapeFunction3(points[i].getX(), points[i].getY()));
            shapeFunctions.set(i,3,  shapeFunction4(points[i].getX(), points[i].getY()));
        }

        trans_matrix = shapeFunctions.transpose();
    }

    public void generateLocalElement(){
        int numOfSidesInEL=4;
        Node [] localNodes = new Node[lpktc*lpktc];

        Point [] integrationPoints = IntegralPoints.getIntegralPoints(lpktc);

        for (int i = 0; i < localNodes.length; i++){
            localNodes[i] = new Node(integrationPoints[i]); }

        int [] Ids= new int[lpktc*lpktc];  Arrays.fill(Ids,0);
        localElement = new Element(Ids, localNodes);
        //2pkt
        double pierwisatek = Math.sqrt(3);
        double zmienna= 1/pierwisatek;
        //3pkt
         double dzielenie=3.0/5.0;
            double pierwiastek3= Math.sqrt(dzielenie);
        //4pkt
        double pierwiastek1= Math.sqrt(6.0/5.0);
            double dziel1=3.0/7.0;
            double dziel2= 2.0/7.0;
            double zmienna1=Math.sqrt(dziel1 -dziel2 * pierwiastek1); //0.33
            double zmienna2=Math.sqrt(dziel1 + dziel2 *pierwiastek1); //0.8

        //1D
        Sides[] sides = new Sides[4];
        if(lpktc==2){
        sides[0] = new Sides( new Node( new Point(-1.0, zmienna)), new Node(new Point(-1.0, -1.0*zmienna)));
        sides[1] = new Sides( new Node( new Point(-1.0*zmienna, -1.0)),  new Node(new Point(zmienna, -1.0)));
        sides[2] = new Sides( new Node( new Point(1.0, -1.0*zmienna)), new Node(new Point(1.0, zmienna)));
        sides[3] = new Sides( new Node( new Point(zmienna, 1.0)),   new Node(new Point(-1.0*zmienna, 1.0))); }
        else if(lpktc==3){
            sides[0] = new Sides( new Node( new Point(-1.0, pierwiastek3)),  new Node(new Point(-1.0, 0)),new Node( new Point(-1.0, -1.0*pierwiastek3))        );
            sides[1] = new Sides( new Node( new Point(-1.0*pierwiastek3, -1.0)), new Node(new Point(0.0, -1.0)),new Node( new Point(pierwiastek3, -1.0))        );
            sides[2] = new Sides( new Node( new Point(1.0, -1.0*pierwiastek3)),   new Node(new Point(1.0, 0.0)), new Node( new Point(1.0, pierwiastek3))      );
            sides[3] = new Sides( new Node( new Point(pierwiastek3, 1.0)),   new Node(new Point(0.0, 1.0)), new Node( new Point(-1.0*pierwiastek3, 1.0))      );
        } else if(lpktc==4){

            sides[0] = new Sides( new Node( new Point(-1.0, zmienna2)),  new Node(new Point(-1., zmienna1)),new Node( new Point(-1.0, -1.0*zmienna1))   ,  new Node( new Point(-1.0, -1.0*zmienna2))   );
            sides[1] = new Sides( new Node( new Point(-1.0*zmienna2, -1.0)), new Node(new Point(-1.0*zmienna1, -1.0)),new Node( new Point(zmienna1, -1.))   ,  new Node( new Point(zmienna2, -1.))   );
            sides[2] = new Sides( new Node( new Point(1.0, -1.0*zmienna2)),   new Node(new Point(1.0, -1.0*zmienna1)), new Node( new Point(1.0, zmienna1))   ,  new Node( new Point(1.0, zmienna2)) );
            sides[3] = new Sides( new Node( new Point(zmienna2, 1.)),   new Node(new Point(zmienna1, 1.0)), new Node( new Point(-1.0*zmienna1, 1.))  ,  new Node( new Point(-1.0*zmienna2, 1.))  );
        }

        localElement.setSides(sides);
        //1D
        for(int i = 0; i < numOfSidesInEL; i++){
            double [][] N_Fun_vals = new double[lpktc][4];
            for(double[] row: N_Fun_vals) { Arrays.fill(row,0.0); }

            for(int j = 0; j < lpktc; j++){
               N_Fun_vals[j][0] =  shapeFunction1(sides[i].getSide()[j].getX(), sides[i].getSide()[j].getY());
                N_Fun_vals[j][1] =  shapeFunction2(sides[i].getSide()[j].getX(), sides[i].getSide()[j].getY());
                N_Fun_vals[j][2] =  shapeFunction3(sides[i].getSide()[j].getX(), sides[i].getSide()[j].getY());
                N_Fun_vals[j][3] = shapeFunction4(sides[i].getSide()[j].getX(), sides[i].getSide()[j].getY());
            }
            sides[i].setShapeFunctionVals(N_Fun_vals);
        }

    }

///////////
public static double shapeFunction1(double ksi, double eta){
    return shapeFunctionBase(ksi, eta, -1., -1.);
    }

    public static double shapeFunction2(double ksi, double eta){
        return shapeFunctionBase(ksi, eta, -1., 1.);
    }

    public static double shapeFunction3(double ksi, double eta){
        return shapeFunctionBase(ksi, eta, 1., 1.);
    }

    public static double shapeFunction4(double ksi, double eta){
        return shapeFunctionBase(ksi, eta, 1., -1.);
    }

    private static double shapeFunctionBase(double ksi, double eta, double etaSign, double ksiSign){
        return 0.25 * ((1 + etaSign * eta) * (1 + ksiSign * ksi));
    }

    //shape functions Eta
    public static double shapeFunction1Eta(double ksi){
        return -0.25 * (1 - ksi);
    }

    public static double shapeFunction2Eta(double ksi) {
        return -0.25 * (1 + ksi);
    }

    public static double shapeFunction3Eta(double ksi){
        return 0.25 * (1 + ksi);
    }

    public static double shapeFunction4Eta(double ksi){
        return 0.25 * (1 - ksi);
    }

    // shape functions Ksi
    public static double shapeFunction1Ksi(double eta){
        return -0.25 * (1 - eta);
    }

    public static double shapeFunction2Ksi( double eta) {
        return 0.25 * (1 - eta);
    }

    public static double shapeFunction3Ksi(double eta){
        return 0.25 * (1 + eta);
    }

    public static double shapeFunction4Ksi(double eta){
        return -0.25 * (1 + eta);
    }



    public Matrix generate_C(double detJ){
        cGlobal = new Matrix(nh, nh);
        double[] waga = IntegralPoints.getWeights2D(lpktc);

        for(int elemIter = 0; elemIter < ne; elemIter++) {

            Element currentElement = (Element)(grid.getEL().get(elemIter));
            cLocal = new Matrix(4,4);

            for (int ipIter = 0; ipIter < lpktc*lpktc; ipIter++) {
                for (int i = 0; i < 4; i++) {
                    for (int j = 0; j < 4; j++) {
                        double cMatrix = c * ro * shapeFunctions.get(ipIter, j) * detJ * trans_matrix.get(i, ipIter) *waga[ipIter] + cLocal.get(i, j); //ok
                        cLocal.set(i, j, cMatrix);
                    }
                }
            }

                for(int i = 0; i < 4; i++){
                    for(int j = 0; j < 4; j++) {
                        int first = currentElement.getIDArray()[i];
                        int second = currentElement.getIDArray()[j];
                        double tempValueC = cGlobal.get(first, second) +  cLocal.get(i,j);
                        cGlobal.set(first, second, tempValueC);
                    }
                }
        }
      return cGlobal;
    }


    public void generate_J_H() {
        hGlobal = new Matrix(nh, nh);
        pGlobal = new Matrix(nh, 1);
        for(int i=0; i<nh;i++){
        pGlobal.set(i, 0, 0.0); }
        Jacobian jacobian;
        double[] dNdx = new double[4];
        double[] dNdy = new double[4];
        double[] wspolrz_X = new double[4];
        double[] wspolrz_Y = new double[4];
        int id; double detJ = 0.;
        double[] node_Temps = new double[4];
        Temp=0.0;


        for(int elemIter = 0; elemIter < ne; elemIter++) {
            Element currentElement = (Element) (grid.getEL().get(elemIter));
            hLocal = new Matrix(4, 4);
            border_cond = new Matrix(4, 4);
            Hbc= new Matrix(4,4);
            p_local= new Matrix(4,1);
            p_local_suma=new Matrix(4,1);
            int numOfNodesInElement=4;


            for (int i = 0; i < numOfNodesInElement; i++) {
                id = currentElement.getIDArray()[i];
                wspolrz_X[i] = ((Node) (grid.getND().get(id))).getX();
                wspolrz_Y[i] = ((Node) (grid.getND().get(id))).getY();
                node_Temps[i] = ((Node) (grid.getND().get(id))).getTemp();
            }

            for (int ipIter = 0; ipIter < lpktc*lpktc; ipIter++){
                double[] waga = IntegralPoints.getWeights2D(lpktc);
                jacobian = new Jacobian(lpktc ,ipIter, wspolrz_X, wspolrz_Y, shapeFunctionsEta, shapeFunctionsKsi);
                  Temp=0.0; int funShape_N=4;

                for (int i = 0; i < funShape_N; i++) {
                    dNdx[i] = jacobian.getJacobianInverseMatrix_final().get(0, 0) * shapeFunctionsKsi.get(ipIter, i) + jacobian.getJacobianInverseMatrix_final().get(0, 1) * shapeFunctionsEta.get(ipIter, i);
                    dNdy[i] = jacobian.getJacobianInverseMatrix_final().get(1, 0) * shapeFunctionsKsi.get(ipIter, i) + jacobian.getJacobianInverseMatrix_final().get(1, 1) * shapeFunctionsEta.get(ipIter, i);
                    Temp += node_Temps[i] * this.shapeFunctions.get(ipIter, i);  }

                detJ = Math.abs(jacobian.getDet());  det_J_global = detJ;

                // MACIERZ  H lokalna
                for (int i = 0; i < 4; i++) {
                    for (int j = 0; j < 4; j++) {
                        double H_matrix_local = hLocal.get(i, j) + k * (dNdx[i] * dNdx[j] + dNdy[i] * dNdy[j])*waga[ipIter] * detJ;
                        hLocal.set(i, j, H_matrix_local);
                    }
                }
            }



            //WARUNKI BRZEGOWE
            for (int bok = 0; bok < currentElement.getNumberOfSideswithBC(); bok++) {

                id = currentElement.getIDOfSideswithBC().get(bok);
                Sides side = currentElement.getSideOfId(id);

                detJ = Math.abs(Math.sqrt(Math.pow((side.getSide()[0].getX() - side.getSide()[1].getX()), 2) + Math.pow((side.getSide()[0].getY() - side.getSide()[1].getY()), 2) ) / 2);
                double[] waga = IntegralPoints.getWeights1D(lpktc);

                for (int i = 0; i < lpktc; i++) {
                    for (int j = 0; j < 4; j++) {
                        for (int k = 0; k < 4; k++) {
                            double border_val = localElement.getSides()[id].getN_Fun_Vals()[i][j]  * localElement.getSides()[id].getN_Fun_Vals()[i][k] *alfa;
                            border_cond.set(j, k, border_val);
                        }

                        double p_local_val= localElement.getSides()[id].getN_Fun_Vals()[i][j] *alfa *waga[i];
                        p_local.set(j,0,p_local_val);
                    }
                    p_local_suma.plusEquals(p_local);

                    border_cond.timesEquals(waga[i]);
                    Hbc.plusEquals(border_cond);
                }

                if(currentElement.getNumberOfSideswithBC() -1 == bok) {
                    p_local_suma.timesEquals(detJ);
                    p_local_suma.timesEquals(-1.0);
                    p_local_suma.timesEquals(t_alfa);
                    Hbc.timesEquals(detJ);
                }
            }


            //agregacja macierzy H i P
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    int first = currentElement.getIDArray()[i];
                    int second = currentElement.getIDArray()[j];
                    double tempValue = hGlobal.get(first, second)  + hLocal.get(i, j) + Hbc.get(i,j);
                    hGlobal.set(first, second, tempValue);
                }
                int first_p = currentElement.getIDArray()[i];
                double temp_p=pGlobal.get(first_p,0) + p_local_suma.get(i,0);
                pGlobal.set(first_p,0,temp_p);
            }
        }
    }



    private GlobalData generate_data(){
        try {
            Gson gson = new Gson();
            GlobalData tempData = gson.fromJson(new FileReader(System.getProperty("user.dir") + "/data/data.txt"), GlobalData.class);

            return tempData;
        } catch (Exception e){
            e.printStackTrace();
        }
        return null;
    }

    public Grid getGrid() {
        return grid;
    }
    public double getDet_J_global() {
        return det_J_global;
    }

    public Matrix getpGlobal() {
        return pGlobal;
    }
    public double getK() {
        return k;
    }

    public void setK(double k) {
        this.k = k;
    }

    public Matrix gethCurrent() {
        return hLocal;
    }
    public double getTau() {
        return tau;
    }
    public void sethCurrent(Matrix hLocal) {
        this.hLocal = hLocal;
    }

    public Matrix gethGlobal() {
        return hGlobal;
    }

    public void sethGlobal(Matrix hGlobal) {
        this.hGlobal = hGlobal;
    }

    public void setShapeFunctionsDerEta(Matrix shapeFunctionsEta) {
        this.shapeFunctionsEta = shapeFunctionsEta;
    }

    public void setShapeFunctionsDerPsi(Matrix shapeFunctionsKsi) {
        this.shapeFunctionsKsi = shapeFunctionsKsi;
    }

    public Matrix getShapeFunctionsEta() {
        return shapeFunctionsEta;
    }

    public Matrix getShapeFunctionsKsi() {
        return shapeFunctionsKsi;
    }

    public double getH() {
        return H;
    }

    public double getT_alfa() {
        return t_alfa;
    }

    public double getDt() {
        return dtau;
    }

    public void setH(double h) {
        H = h;
    }

    public double getB() {
        return B;
    }

    public void setB(double b) {
        B = b;
    }

    public int getnH() {
        return nH;
    }

    public void setnH(int nH) {
        this.nH = nH;
    }

    public int getnB() {
        return nB;
    }

    public void setnB(int nB) {
        this.nB = nB;
    }

    public int getNh() {
        return nh;
    }

    public void setNh(int nh) {
        this.nh = nh;
    }

    public int getNe() {
        return ne;
    }

    public void setNe(int ne) {
        this.ne = ne;
    }

    public double getDy() {
        return dy;
    }

    private void setDy(double dy) {
        this.dy = dy;
    }

    public double getDx() {
        return dx;
    }

    private void setDx(double dx) {
        this.dx = dx;
    }

    public int getLpktc() {
        return lpktc;
    }

    public double getC() {
        return c;
    }

    public double getRo() {
        return ro;
    }

    public double getT0() {
        return t0;
    }

    public double getAlfa() {
        return alfa;
    }

    public double getDtau() {
        return dtau;
    }
}
