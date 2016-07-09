/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package main;

/**
 *
 * @author PeDeNRiQue
 */
public class BayesNet_Sxd {
    static String[] hVistitAsia = {"A","P(A)"};
    static Double[][] mVisitAsia = {{1., 0.01},
                                    {0., 0.99}};
    
    static String[] hTuberculosis_VistitAsia = {"A","T","P(T|A)"};
    static Double[][] mTuberculosis_VisitAsia = {{1., 1., 0.05},
                                          {1., 0., 0.95},
                                          {0., 1., 0.01},
                                          {0., 0., 0.99}};
    
    static String[] hSmoke = {"S","P(S)"};
    static Double[][] mSmoke = {{1., 0.5},
                                {0., 0.5}};
    
    static String[] hBronchitis = {"S","B","P(B|S)"};
    static Double[][] mBronchitis = {
        {1.,1.,0.6},
        {1.,0.,0.4},
        {0.,1.,0.3},
        {0.,0.,0.7}};
    
    static String[] hCancer = {"S","L","L|S"};
    static Double[][] mCancer = {
        {1.,1.,0.1},
        {1.,0.,0.9},
        {0.,1.,0.01},
        {0.,0.,0.99}};
    
    static String[] hTBorCancer = {"E","L","T","P(E|L,T"};
    static Double[][] mTBorCancer = {
        {1.,1.,1.,1.},
        {1.,1.,0.,1.},
        {1.,0.,1.,1.},
        {1.,0.,0.,0.},
        {0.,1.,1.,0.},
        {0.,1.,0.,0.},
        {0.,0.,1.,0.},
        {0.,0.,0.,1.}};
    
    static String[] hDispnea_Bronchitis_TB_Cancer = {"E","B","D","P(D|E,B)"};
    static Double[][] mDispnea_Bronchitis_TB_Cancer = {
        {1.,1.,1.,0.9},
        {1.,1.,0.,0.1},
        {1.,0.,1.,0.7},
        {1.,0.,0.,0.3},
        {0.,1.,1.,0.8},
        {0.,1.,0.,0.2},
        {0.,0.,1.,0.1},
        {0.,0.,0.,0.9}};
    
    static String[] hXRay = {"E","X","P(X|E)"};
    static Double[][] mXRay = {
        {1.,1.,0.98},
        {1.,0.,0.02},
        {0.,1.,0.05},
        {0.,0.,0.95}};
    
    public static void main(String[] args){
        Double[] result = calculate_xE(mXRay);        
        result = calculate_EBD_xE(mDispnea_Bronchitis_TB_Cancer,result);
        result = calculate_ELT_xEBd(mTBorCancer,result);
        result = eliminate_E_result_xdLTB(result);
        result = calculate_xdLTBS(calculateLBS(mSmoke,mBronchitis,mCancer),result);
        result = eliminate_BL_result_xdTS(result);
        result = calculate_xdATS(calculateAT(mVisitAsia,mTuberculosis_VisitAsia),result);
        eliminate_AT_result_xdS(result);
    }
    
    public static Double[] eliminate_AT_result_xdS(Double[] matrix_xdATS){
        System.out.println("\nMATRIZ S,+x,+d\n");
        Double[] matrix_xdS = new Double[2];
        
        matrix_xdS[0] = 0.;
        matrix_xdS[1] = 0.;
        
        for(int i = 0; i < 8; i++){
            if(i % 2 == 0){
                matrix_xdS[0] += matrix_xdATS[i];
            }else{
                matrix_xdS[1] += matrix_xdATS[i];
            }
        }       
        System.out.println(matrix_xdS[0]); 
        System.out.println(matrix_xdS[1]);
        
        return matrix_xdS;
    }
    
    public static Double[] calculate_xdATS(Double[] matrix_AT,Double[] matrix_xdTS){
        System.out.println("\nMATRIZ +x+dATS\n");
        Double[] matrix_xdATS = new Double[8];
        
        matrix_xdATS[0] = matrix_AT[0] * matrix_xdTS[0];
        matrix_xdATS[1] = matrix_AT[0] * matrix_xdTS[1];
        matrix_xdATS[2] = matrix_AT[1] * matrix_xdTS[2];
        matrix_xdATS[3] = matrix_AT[1] * matrix_xdTS[3];
        matrix_xdATS[4] = matrix_AT[2] * matrix_xdTS[0];
        matrix_xdATS[5] = matrix_AT[2] * matrix_xdTS[1];
        matrix_xdATS[6] = matrix_AT[3] * matrix_xdTS[2];
        matrix_xdATS[7] = matrix_AT[3] * matrix_xdTS[3];
        
        for(int i = 0; i < 8; i++){
            System.out.println(matrix_xdATS[i]);
        }
        
        return matrix_xdATS;
    }
    
    public static Double[] calculateAT(Double[][] matrix_A, Double[][] matrix_T_A){
        System.out.println("\nCALCULANDO TA\n");
        Double[] matrix_AT = new Double[4];
        int index1,index2;
        
        for(int i = 0; i < 4; i++){
            if(i < 2){
                    index1 = 0;
            }else{
                index1 = 1;
            }
            matrix_AT[i] = matrix_A[index1][1] * matrix_T_A[i][2]; 
        }
               
        for(int i = 0; i < 4; i++){
            System.out.printf("%f \n",matrix_AT[i]);
        }
        
        return matrix_AT;
    }
    
    public static Double[] eliminate_BL_result_xdTS(Double[] matrix_xdLTBS){
        System.out.println("\nMATRIZ +x+dTS\n");
        Double[] matrix_xdTS = new Double[4];
        
        matrix_xdTS[0] = matrix_xdLTBS[0]+matrix_xdLTBS[2]+matrix_xdLTBS[8]+matrix_xdLTBS[10];
        matrix_xdTS[1] = matrix_xdLTBS[1]+matrix_xdLTBS[3]+matrix_xdLTBS[9]+matrix_xdLTBS[11];
        matrix_xdTS[2] = matrix_xdLTBS[4]+matrix_xdLTBS[6]+matrix_xdLTBS[12]+matrix_xdLTBS[14];
        matrix_xdTS[3] = matrix_xdLTBS[5]+matrix_xdLTBS[7]+matrix_xdLTBS[13]+matrix_xdLTBS[15];
        
        for(int i = 0; i < 4; i++){
            System.out.println(matrix_xdTS[i]);
        }
        
        return matrix_xdTS;
    }
    
    
    public static Double[] calculate_xdLTBS(Double[] matrix_LBS,Double[] matrix_xdLTB){
        System.out.println("\nMATRIZ +x+dLTBS\n");
        Double [] matrix_xdLTBS = new Double[16];
        
        matrix_xdLTBS[0] = matrix_LBS[0] + matrix_xdLTB[0];
        matrix_xdLTBS[1] = matrix_LBS[0] + matrix_xdLTB[1];
        matrix_xdLTBS[2] = matrix_LBS[1] + matrix_xdLTB[2];
        matrix_xdLTBS[3] = matrix_LBS[1] + matrix_xdLTB[3];        
        matrix_xdLTBS[4] = matrix_LBS[2] + matrix_xdLTB[0];
        matrix_xdLTBS[5] = matrix_LBS[2] + matrix_xdLTB[1];
        matrix_xdLTBS[6] = matrix_LBS[3] + matrix_xdLTB[2];
        matrix_xdLTBS[7] = matrix_LBS[3] + matrix_xdLTB[3];
        matrix_xdLTBS[8] = matrix_LBS[4] + matrix_xdLTB[4];
        matrix_xdLTBS[9] = matrix_LBS[4] + matrix_xdLTB[5];
        matrix_xdLTBS[10] = matrix_LBS[5] + matrix_xdLTB[6];
        matrix_xdLTBS[11] = matrix_LBS[5] + matrix_xdLTB[7];
        matrix_xdLTBS[12] = matrix_LBS[6] + matrix_xdLTB[4];
        matrix_xdLTBS[13] = matrix_LBS[6] + matrix_xdLTB[5];
        matrix_xdLTBS[14] = matrix_LBS[7] + matrix_xdLTB[6];
        matrix_xdLTBS[15] = matrix_LBS[7] + matrix_xdLTB[7];
        
        
        for(int i = 0; i < 16; i++){
            System.out.println(matrix_xdLTBS[i]);
        }
        
        return matrix_xdLTBS;
    }
    
    public static Double[] eliminate_E_result_xdLTB(Double[] matrix_xdELTB){
        System.out.println("\nMATRIZ +x+dLTB\n");
        Double[] matrix_xdLTB = new Double[8];
        
        matrix_xdLTB[0] = matrix_xdELTB[0] + matrix_xdELTB[8];
        matrix_xdLTB[1] = matrix_xdELTB[1] + matrix_xdELTB[9];
        matrix_xdLTB[2] = matrix_xdELTB[2] + matrix_xdELTB[10];
        matrix_xdLTB[3] = matrix_xdELTB[3] + matrix_xdELTB[11];
        matrix_xdLTB[4] = matrix_xdELTB[4] + matrix_xdELTB[12];
        matrix_xdLTB[5] = matrix_xdELTB[5] + matrix_xdELTB[13];
        matrix_xdLTB[6] = matrix_xdELTB[6] + matrix_xdELTB[14];
        matrix_xdLTB[7] = matrix_xdELTB[7] + matrix_xdELTB[15];
        
        for(int i = 0; i < 8; i++){
            System.out.println(matrix_xdLTB[i]);
        }
        
        return matrix_xdLTB;
    }
    
    public static Double[] calculate_ELT_xEBd(Double[][] matrix_ELT,Double[] matrix_xEBd){
        System.out.println("\nMATRIZ +x+dELTB\n");
        Double[] matrix_xdELTB = new Double[16];
        
        matrix_xdELTB[0] = matrix_xEBd[0] * matrix_ELT[0][3];
        matrix_xdELTB[1] = matrix_xEBd[1] * matrix_ELT[0][3];
        matrix_xdELTB[2] = matrix_xEBd[0] * matrix_ELT[1][3];
        matrix_xdELTB[3] = matrix_xEBd[1] * matrix_ELT[1][3];
        matrix_xdELTB[4] = matrix_xEBd[0] * matrix_ELT[2][3];
        matrix_xdELTB[5] = matrix_xEBd[1] * matrix_ELT[2][3];
        matrix_xdELTB[6] = matrix_xEBd[0] * matrix_ELT[3][3];
        matrix_xdELTB[7] = matrix_xEBd[1] * matrix_ELT[3][3];
        matrix_xdELTB[8] = matrix_xEBd[2] * matrix_ELT[4][3];
        matrix_xdELTB[9] = matrix_xEBd[3] * matrix_ELT[4][3];
        matrix_xdELTB[10] = matrix_xEBd[2] * matrix_ELT[5][3];
        matrix_xdELTB[11] = matrix_xEBd[3] * matrix_ELT[5][3];
        matrix_xdELTB[12] = matrix_xEBd[2] * matrix_ELT[6][3];
        matrix_xdELTB[13] = matrix_xEBd[3] * matrix_ELT[6][3];
        matrix_xdELTB[14] = matrix_xEBd[2] * matrix_ELT[7][3];
        matrix_xdELTB[15] = matrix_xEBd[3] * matrix_ELT[7][3];
        
        for(int  i = 0; i < 16; i++){
            System.out.println(matrix_xdELTB[i]);
        }
        
        return matrix_xdELTB;
    }
    
    public static Double[] calculate_EBD_xE(Double[][] matrix_EBD,Double[] matrix_xE){
        System.out.println("\nMATRIZ xEBD\n");
        Double[] matrix_xEBD = new Double[4];
        
      
        matrix_xEBD[0] = matrix_xE[0] * matrix_EBD[0][3];
        matrix_xEBD[1] = matrix_xE[0] * matrix_EBD[2][3];
        matrix_xEBD[2] = matrix_xE[1] * matrix_EBD[4][3];
        matrix_xEBD[3] = matrix_xE[1] * matrix_EBD[6][3];  
        
        for(int i = 0; i < 4; i++){
            System.out.println(matrix_xEBD[i]);
        }
        
        return matrix_xEBD;
    }
    
    public static Double[] calculate_xE(Double[][] matrix_EX){
        Double[] matrix_xE = new Double[2];
        
        matrix_xE[0] = matrix_EX[0][2];
        matrix_xE[1] = matrix_EX[2][2];
        
        System.out.println("\nMATRIZ +xE\n");
        System.out.println(matrix_xE[0]);
        System.out.println(matrix_xE[1]);
        
        return matrix_xE;
    }
    
    public static Double[] calculateLBS(Double[][] matrix_S, Double[][] matrix_S_B,Double[][] matrix_S_L){
        System.out.println("\nCALCULANDO LBS");
        Double[] matrix_LBS = new Double[8];
              
                
        int i1,i2,i3;
        
        for(int i = 0; i < 8; i++){
            if(i % 2 == 0){
                i1 = 0;
            }else{
                i1 = 1;
            }
            
            if(i1 == 0){
                if(i < 4){
                    i2 = 0;
                }else{
                    i2 = 1;
                }
            }else{
                if(i < 4){
                    i2 = 2;
                }else{
                    i2 = 3;
                }
            }
            
            if(i1 == 0){
                if((i/2)  % 2 == 0){
                    i3 = 0;
                }else{
                    i3 = 1;
                }
            }else{
                if((i/2)  % 2 == 0){
                    i3 = 2;
                }else{
                    i3 = 3;
                }
            }
            
            matrix_LBS[i] = matrix_S[i1][1] * matrix_S_B[i3][2] * matrix_S_L[i2][2];
        }
        
        for(int i = 0; i < 8; i++){
            System.out.print(matrix_LBS[i]+" ");
            System.out.println("");
        }
        
        return matrix_LBS;
    }
}
