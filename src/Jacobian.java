import Jama.Matrix;


public class Jacobian {
    private Matrix J;
    private Matrix jInverted;
    private Matrix shapeFunctionsEta;
    private Matrix shapeFunctionsKsi;
    private double det;

    public Jacobian( int lpktc ,int integrationPoint, double x[], double y[], Matrix shapeFunctionsEta, Matrix shapeFunctionsKsi){
        this.shapeFunctionsEta = shapeFunctionsEta;
        this.shapeFunctionsKsi = shapeFunctionsKsi;
        J = new Matrix(2,2);

         double dXdKsi=0.0, dYdKsi=0.0, dXdEta=0.0, dYdEta=0.0;

             dXdKsi = this.shapeFunctionsKsi.get(integrationPoint, 0) * x[0] + this.shapeFunctionsKsi.get(integrationPoint, 1) * x[1]
                    + this.shapeFunctionsKsi.get(integrationPoint, 2) * x[2] + this.shapeFunctionsKsi.get(integrationPoint, 3) * x[3];

             dYdKsi = this.shapeFunctionsKsi.get(integrationPoint, 0) * y[0] + this.shapeFunctionsKsi.get(integrationPoint, 1) * y[1]
                    + this.shapeFunctionsKsi.get(integrationPoint, 2) * y[2] + this.shapeFunctionsKsi.get(integrationPoint, 3) * y[3];

             dXdEta = this.shapeFunctionsEta.get(integrationPoint, 0) * x[0] + this.shapeFunctionsEta.get(integrationPoint, 1) * x[1]
                    + this.shapeFunctionsEta.get(integrationPoint, 2) * x[2] + this.shapeFunctionsEta.get(integrationPoint, 3) * x[3];

            dYdEta = this.shapeFunctionsEta.get(integrationPoint, 0) * y[0] + this.shapeFunctionsEta.get(integrationPoint, 1) * y[1]
                    + this.shapeFunctionsEta.get(integrationPoint, 2) * y[2] + this.shapeFunctionsEta.get(integrationPoint, 3) * y[3];


        J.set(0,0,dXdKsi);
        J.set(0,1,dYdKsi);
        J.set(1,0,dXdEta);
        J.set(1,1,dYdEta);


        det = J.det();
        jInverted = J.transpose();
        jInverted.set(0,1, jInverted.get(0,1)*(-1.0) );
        jInverted.set(1,0, jInverted.get(1,0)*(-1.0) );
    }

    public Matrix getJacobianInverseMatrix_final(){
        return jInverted.times(1.0/det);
    }

    public double getDet() {
        return det;
    }
}
