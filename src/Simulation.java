import Jama.Matrix;

import java.io.PrintWriter;
import java.util.Arrays;
import  java.io.PrintWriter.*;
public class Simulation {

    public static void main(String[] args) {

        GlobalData x = new GlobalData(true);
        double dTAU = x.getDt();
        double Tau=x.getTau();
        int it=1;
        double it_div=Tau/dTAU;
        double[] help_array=new double[(int)it_div];
        double[] min_temp_array=  new double[(int)it_div];
        double[] max_temp_array=  new double[(int)it_div] ;
        double t0 = x.getT0();
        int l_nodow = x.getNh();
        double[] t_table= new double[l_nodow];
        Arrays.fill(t_table,t0);
        double div_Tau = 1.0 / dTAU;


        for(int time=0; time<Tau ; time+=dTAU ){

            x.generate_J_H();
            double detJ = x.getDet_J_global();

            Matrix C_global = new Matrix(l_nodow, l_nodow);
            Matrix H_global = new Matrix(l_nodow, l_nodow);
            Matrix C_div_dt = new Matrix(l_nodow, l_nodow);
            double[] P_Array = new double[l_nodow];
            Arrays.fill(P_Array,0);
            Matrix vector_c_t = new Matrix(l_nodow, 1);

            Matrix P_Global = x.getpGlobal();
            C_global = x.generate_C(detJ);
            H_global = x.gethGlobal();

            // Martix [H] = [H]+[C]/dT
            C_global.timesEquals(div_Tau);
            H_global.plusEquals(C_global); //H dash

            // and {P} = {P}+{[C]/dT}*{T0}
            C_div_dt = C_global;         //= [C]/dT
            for (int j = 0; j < l_nodow; j++) {
                for (int k = 0; k < l_nodow; k++) {
                    P_Array[j] += C_div_dt.get(j, k)*t_table[k]; //{P} = {P}+ ....*{T0}
                }
            }

            for (int j = 0; j < l_nodow; j++) {
                vector_c_t.set(j, 0, P_Array[j]);
            }

            P_Global.timesEquals(-1);
            P_Global.plusEquals(vector_c_t);  //P dash
            

            t_table=generate_solution(l_nodow,P_Global,H_global);

            for(int i = 0; i < x.getNh(); i++){
                ((Node)(x.getGrid().getND().get(i))).setTemp(t_table[i]);
            }

                help_array=t_table.clone();
                Arrays.sort(help_array);
                min_temp_array[ it-1 ]= help_array[0];
                max_temp_array[it -1] = help_array[l_nodow-1];
            it++;
        }

        for(int i=0; i<20;i++){
            System.out.println(i +" Min: " + min_temp_array[i] + " Max "+ max_temp_array[i]);
        }

    }

    public static double[] generate_solution(int l_nodow, Matrix P_Global, Matrix H_global){
       double max=0;
       int wiersz=0;
       int wiersz2=0;
        double[][] Hg_arr= new double[l_nodow][l_nodow];
        double[] Pg_arr= new double[l_nodow];
        double[] Result= new double[l_nodow];

        //dane z Matrix H i wektora P do tablic
        for(int i=0; i<l_nodow; i++){
            for(int j=0; j<l_nodow; j++){
                Hg_arr[i][j]=H_global.get(i,j);
            }
            Pg_arr[i]=P_Global.get(i,0);
        }

//szukanie najwiekszego elementu i zamiana wierszy
       for(int i=0; i<l_nodow; i++){
           if (Hg_arr[i][i] == 0)
           {
               wiersz=i;
               for (int j = 0; j < l_nodow; j++)
               {
                   if (Math.abs(Hg_arr[j][i]) > max)
                   {
                       max = Hg_arr[j][i];
                       wiersz2=j;
                   }
               }
           }

           if (max != 0)
           {
               for (int k = 0; k < l_nodow; k++)
               {    //swap
                   double temp1= Hg_arr[wiersz][k];
                   double temp2=Hg_arr[wiersz2][k];
                   Hg_arr[wiersz][k]=temp2;
                   Hg_arr[wiersz2][k]=temp1;
               }
               //swap
               double val1= Pg_arr[wiersz];
               double val2=Pg_arr[wiersz2];
               Pg_arr[wiersz]=val2;
               Pg_arr[wiersz2]=val1;
               max=0;
               //System.out.println("Wystapilo 0 na przekatnej, macierz po zamianie wierszy: ");
           }

       }

          // Operacje przeksztalcania macierzy metoda Gaussa -eliminacja wpolczynnikow
           for (int k = 0; k < l_nodow - 1; k ++)
           {
               for (int i = k + 1; i < l_nodow; i++)
               {
                   Hg_arr[i][k] /= Hg_arr[k][k];
                   for (int j = k + 1; j < l_nodow + 1; j++)
                   {
                       if (j == l_nodow)
                           Pg_arr[i] -= Hg_arr[i][k] * Pg_arr[k];
                       else
                           Hg_arr[i][j] -= Hg_arr[i][k] * Hg_arr[k][j];
                   }
                   Hg_arr[i][k] = 0;
               }
           }

           //Obliczanie wartosci niewiadomych -> Result
           for ( int i = l_nodow - 1; i >= 0; i --)
           {
               Result[i] = Pg_arr[i];
               for (int j = i + 1; j < l_nodow; j++)
               {
                   if (j == l_nodow)
                       Result[i] -= Pg_arr[i] * Result[j];
                   else
                       Result[i] -= Hg_arr[i][j] * Result[j];
               }
               Result[i] /= Hg_arr[i][i];
           }
           return Result;
    }
}