public class IntegralPoints {

    private int lpktc;
    private double [][] tab2p = null;
    private double [][] tab3p = null;
    private double [][] tab4p = null;

    public IntegralPoints(int lpktc) {
    if(lpktc==2) {
        tab2p = new double[2][2];
        tab2p[0][0] = -1 / Math.sqrt(3);
        tab2p[0][1] = 1.;
        tab2p[1][0] = 1 / Math.sqrt(3);
        tab2p[1][1] = 1.;
    }else
        if(lpktc==3) {
        tab3p = new double[3][2];

        tab3p[0][0] = -0.7745;  tab3p[0][1] = 5 / (double) 9;
        tab3p[1][0] = 0.;       tab3p[1][1] = 8 / (double) 9;
        tab3p[2][0] = 0.7745;   tab3p[2][1] = 5 / (double) 9;
    }else
        if (lpktc==4){
        tab4p = new double[4][2]; //ok
        tab4p[0][0]=-0.861136;  tab4p[0][1]=0.347855 ;
        tab4p[1][0]=-0.339981; tab4p[1][1]=0.652145;
        tab4p[2][0]=0.339981; tab4p[2][1]=0.652145;
        tab4p[3][0]=0.861136;  tab4p[3][1]=0.347855 ; }
    }

    public static Point[] getIntegralPoints(int lpktc ){
        Point [] points = new Point[lpktc*lpktc];
        if( lpktc ==2)
        {
            double pierwiastek = Math.sqrt(3.0);
            double zmienna = 1.0/pierwiastek;
             //ksi = X, eta =Y
             points[0] = new Point(-1.0*zmienna, -1.0*zmienna);
             points[1] = new Point(zmienna, -1.0*zmienna);
             points[2] = new Point(zmienna, zmienna);
             points[3] = new Point(-1.0*zmienna, zmienna);
            return points;
        }
        else if(lpktc==3){
            double dzielenie=3.0/5.0;
            double pierwiastek= Math.sqrt(dzielenie);
            points[0] = new Point(-1.0*pierwiastek, -1.0*pierwiastek);
            points[1] = new Point(0.0, -1.0*pierwiastek);
            points[2] = new Point(pierwiastek, -1.0*pierwiastek);

            points[3] = new Point(-1.0*pierwiastek, 0.0);
            points[4] = new Point(0.0, 0.0);
            points[5] = new Point(pierwiastek, 0.0);

            points[6] = new Point(-1.0*pierwiastek, pierwiastek);
            points[7] = new Point(0.0, pierwiastek);
            points[8] = new Point(pierwiastek, pierwiastek);
            return points;
        }else if(lpktc==4){
            double pierwiastek1= Math.sqrt(6.0/5.0);
            double dziel1=3.0/7.0;
            double dziel2= 2.0/7.0;
            double zmienna=Math.sqrt(dziel1 -dziel2 * pierwiastek1);
            double zmienna2=Math.sqrt(dziel1 + dziel2 *pierwiastek1);

            points[0]= new Point(-1.0*zmienna2,-1.0*zmienna2);
            points[1]= new Point(-1.0*zmienna,-1.0*zmienna2);
            points[2]= new Point(zmienna,-1.0*zmienna2);
            points[3]= new Point(zmienna2,-1.0*zmienna2);

            points[4]= new Point(-1.0*zmienna2,-1.0*zmienna);
            points[5]= new Point(-1.0*zmienna,-1.0*zmienna);
            points[6]= new Point(zmienna,-1.0*zmienna);
            points[7]= new Point(zmienna2,-1.0*zmienna);

            points[8]= new Point(-1.0*zmienna2,zmienna);
            points[9]= new Point(-1.0*zmienna,zmienna);
            points[10]= new Point(zmienna,zmienna);
            points[11]= new Point(zmienna2,zmienna);

            points[12]= new Point(-1.0*zmienna2,zmienna2);
            points[13]= new Point(-1.0*zmienna,zmienna2);
            points[14]= new Point(zmienna,zmienna2);
            points[15]= new Point(zmienna2,zmienna2);
            return points;
        }

        return points;
    }

    public static double[] getWeights1D(int lpktc){

        double [] waga= new double[lpktc];
        if( lpktc ==2)
        {
            waga[0]=1.0;
            waga[1]=1.0;

        }else if(lpktc==3){
            waga[0]=(5/(double)9);
            waga[1]=(8/(double)9);
            waga[2]=(5/(double)9);
        }
        else if(lpktc==4){
            double pierwiastek= Math.sqrt(30.0);
            double zmienna= (18.0+ pierwiastek)/36.0; //0.65
            double zmienna2=(18.0- pierwiastek)/36.0; //0.347
             waga[0]=zmienna2;
             waga[1]=zmienna;
             waga[2]=zmienna;
             waga[3]=zmienna2;
        }
            return waga;
    }



    public static double[] getWeights2D(int lpktc){

        double [] waga= new double[lpktc*lpktc];
        if( lpktc ==2)
        {
            waga[0]=1.0;
            waga[1]=1.0;
            waga[2]=1.0;
            waga[3]=1.0;

        }else if(lpktc==3){
            waga[0]= (5/(double)9) * (5/(double)9) ;
            waga[1]= (8/(double)9) * (  5/(double)9 );
            waga[2]= (5/(double)9)*( 5/(double)9 );

            waga[3]=(5/(double)9) * (8/(double)9)  ;
            waga[4]= (8/(double)9) * (8/(double)9) ;
            waga[5]=(5/(double)9) * (8/(double)9) ;

            waga[6]= (5/(double)9) * (5/(double)9) ;
            waga[7]= (8/(double)9) * (  5/(double)9 );
            waga[8]= (5/(double)9)*( 5/(double)9 );

        }
        else if(lpktc==4){
            double pierwiastek= Math.sqrt(30.0);
            double zmienna= (18.0+ pierwiastek)/36.0; //0.65
            double zmienna2=(18.0- pierwiastek)/36.0; //0.347

            waga[0]=zmienna2*zmienna;
            waga[1]=zmienna *zmienna2;
            waga[2]=zmienna*zmienna2;
            waga[3]=zmienna2 *zmienna2;

            waga[4]=zmienna2* zmienna;
            waga[5]=zmienna* zmienna;
            waga[6]=zmienna* zmienna;
            waga[7]=zmienna2 * zmienna;

            waga[8]=zmienna2 * zmienna;
            waga[9]=zmienna* zmienna;
            waga[10]=zmienna* zmienna;
            waga[3]=zmienna2 * zmienna;

            waga[11]=zmienna2 *zmienna2;
            waga[12]=zmienna*zmienna2;
            waga[13]=zmienna*zmienna2;
            waga[14]=zmienna2 *zmienna2;
        }
        return waga;
    }

    public double[][] getTab2p() {
        return tab2p;
    }

    public double[][] getTab3p() {
        return tab3p;
    }

    public double[][] getTab4p(){
        return tab4p;
    }
}
