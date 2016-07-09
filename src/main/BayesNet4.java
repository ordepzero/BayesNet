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
public class BayesNet4 {
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
    
    public static void main(String[] args) {
        
        Double[] result = calculateTA(mVisitAsia,mTuberculosis_VisitAsia);
        result = eliminateA_resultT(result);
        result = normalize(result);
        result = calculateELT(mTBorCancer,result);
        result = eliminateT_resultEL(result);
        result = normalize(result);
        result = calculate_xEL(mXRay,result);
        result = calculate_XELDB(mDispnea_Bronchitis_TB_Cancer,result);
        result = eliminateE_evidenceD_resultXLdB(result);
        result = normalize(result);
        result = calculateXLBS(calculateLBS(mSmoke,mBronchitis,mCancer),result);
        result = normalize(result);
        calculateSxd(result);
    }
    
    public static Double[] normalize(Double[] matrix){
        Double total = 0.;
        
        for(int i = 0; i < matrix.length; i++){
            total += matrix[i];
        }
        for(int i = 0; i < matrix.length; i++){
            
            matrix[i] = matrix[i]/total;
        }
        
        return matrix;
    }
    
    public static void calculateSxd(Double[] matrix_xdLBS){
        //System.out.println("\nCALCULANDO P(S,+x,+d)\n");
        Double[] matrix_Sxd = new Double[2];
        
        matrix_Sxd[0] = 0.;
        matrix_Sxd[1] = 0.;
        
        for(int i = 0; i < 8; i++){
            if(i % 2 == 0){
                matrix_Sxd[0] += matrix_xdLBS[i];
            }else{
                matrix_Sxd[1] += matrix_xdLBS[i];
            }
        }
        
        
        
        System.out.println("S | P(S,+x,+d)");
        System.out.println("+s| "+matrix_Sxd[0]);
        System.out.println("-s| "+matrix_Sxd[1]);
    }
    
    public static Double[] calculateXLBS(Double[] matrix_LBS,Double[] matrix_XLdB){
        //System.out.println("\nCALCULANDO P(X,L,d,B,S)\n");
        Double[] matrix_XLBS = new Double[8];
        
        for(int i = 0; i < 8; i++){
            matrix_XLBS[i] = matrix_LBS[i] * matrix_XLdB[i / 2];
            //System.out.printf("%f\n",matrix_XLBS[i]);
        }
        
        return matrix_XLBS;
    }
    
    public static Double[] calculateLBS(Double[][] matrix_S, Double[][] matrix_S_B,Double[][] matrix_S_L){
        //System.out.println("\nCALCULANDO LBS");
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
        
//        for(int i = 0; i < 8; i++){
//           System.out.printf("%f\n",matrix_LBS[i]);
//        }
        
        return matrix_LBS;
    }
    
    public static Double[] eliminateE_evidenceD_resultXLdB(Double[] matrix_XELDB){
        //System.out.println("\nELIMINANDO E E EVIDENCIA D -> P(X,L,d,B)\n");
        Double[] matrix_XLdB = new Double[4];
        
        matrix_XLdB[0] = matrix_XELDB[0];
        matrix_XLdB[1] = matrix_XELDB[1];
        matrix_XLdB[2] = matrix_XELDB[4];
        matrix_XLdB[3] = matrix_XELDB[5];
        
//        for(int i = 0; i < 4; i++){
//            System.out.printf("%f\n",matrix_XLdB[i]);
//        }
        
        return matrix_XLdB;
    }    
    
    public static Double[] calculate_XELDB(Double[][] matrix_DEB,Double[] matrix_xEL){
        //System.out.println("\nCALCULANDO P(X,E,L,D,B)\n");
        Double[]  matrix_XELDB = new Double[8];
        int i1 = 0;
        
        for(int i = 0; i < 8; i++){
            if(i < 4){
                
                if(i == 0){
                    i1 = 0;
                }else if(i == 1){
                    i1 = 2;
                }else if(i == 2){
                    i1 = 1;
                }else if(i == 3){
                    i1 = 3;
                }
                matrix_XELDB[i] = matrix_xEL[0] * matrix_DEB[i1][3];
                    
            }else{
                if(i - 4 == 0){
                    i1 = 0;
                }else if(i - 4 == 1){
                    i1 = 2;
                }else if(i - 4 == 2){
                    i1 = 1;
                }else if(i - 4 == 3){
                    i1 = 3;
                }
                matrix_XELDB[i] = matrix_xEL[1] * matrix_DEB[i1][3];
            }
        }
        
//        for(int i = 0; i < 8; i++){
//            System.out.printf("%f \n",matrix_XELDB[i]);
//        }
        
        return matrix_XELDB;
    }
    
    public static Double[] calculate_xEL(Double[][] matrix_EX,Double[] matrixEL){
        //System.out.println("\nCALCULANDO P(X,E,L)\n");
        Double[] matrix_XEL = new Double[2];
        
        for(int i = 0; i < 2; i++){
            if(i % 2 == 0){
                matrix_XEL[i] = matrix_EX[0][2] * matrixEL[0];
            }else{
                matrix_XEL[i] = matrix_EX[0][2] * matrixEL[1];
            }
        }
        
//        for(int  i = 0; i < 2; i++){
//            System.out.printf("%f\n",matrix_XEL[i]);
//        }
            
        return matrix_XEL;
    }
    
    public static Double[] eliminateT_resultEL(Double[] matrix_ELT){
        //System.out.println("\nELIMINANDO T -> P(E,L)\n");
        Double[] matrix_EL = new Double[2];
        
        matrix_EL[0] = matrix_ELT[0] + matrix_ELT[1];
        matrix_EL[1] = matrix_ELT[2] + matrix_ELT[3];
        
//        for(int i = 0; i < 2; i++){
//            System.out.printf("%f\n",matrix_EL[i]);
//        }
        
        return matrix_EL;
    }
    
    public static Double[] calculateELT(Double[][] matrixE_LT, Double[] matrix_T){
        //System.out.println("\nCALCULANDO P(E,L,T)\n");
        Double[] matrix_ELT = new Double[4];
        
        for(int i = 0; i < 4; i++){
            if(i % 2 == 0){
                matrix_ELT[i] = matrixE_LT[i][3] * matrix_T[0];
            }else{
                matrix_ELT[i] = matrixE_LT[i][3] * matrix_T[1];
            }
        }
//        for(int i = 0; i < 4; i++){
//            System.out.printf("%f\n",matrix_ELT[i]);
//        }
        return matrix_ELT;
    }
    
    public static Double[] eliminateA_resultT(Double[] matrix_TA){
        //System.out.println("\nELIMINANDO A -> P(T)");
        Double[] matrix_T = new Double[2];
        matrix_T[0] = 0.;
        matrix_T[1] = 0.;
        
        for(int i = 0; i < 4; i++){
            if(i % 2 == 0){
                matrix_T[0] += matrix_TA[i];
            }else{
                matrix_T[1] += matrix_TA[i];
            }
        }
            
//        for(int i = 0; i < 2; i++){
//            System.out.printf("%f\n",matrix_T[i]);
//        }
        
        return matrix_T;
    }
    
    public static Double[] calculateTA(Double[][] matrix_A, Double[][] matrix_T_A){
        //System.out.println("\nCALCULANDO TA\n");
        Double[] matrix_TA = new Double[4];
        int index1,index2;
        
        for(int i = 0; i < 4; i++){
            if(i < 2){
                    index1 = 0;
            }else{
                index1 = 1;
            }
            matrix_TA[i] = matrix_A[index1][1] * matrix_T_A[i][2]; 
        }
               
//        for(int i = 0; i < 4; i++){
//            System.out.printf("%f \n",matrix_TA[i]);
//        }
        
        return matrix_TA;
    }
}
