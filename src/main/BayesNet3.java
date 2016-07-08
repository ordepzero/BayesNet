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
public class BayesNet3 {
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
        {1.,0.,0.,0.}};
    
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
        
        Double[] result = calculateLBS(mSmoke,mBronchitis,mCancer);
        result = evidenceS_resultLBs(result);
        result = calculateASTLB(calculateTA(mVisitAsia,mTuberculosis_VisitAsia),result);
        result = calculateASETLB(mTBorCancer,result);
        result = eliminateTL_resultASEB(result);
        result = calculateASEBD(mDispnea_Bronchitis_TB_Cancer,result);
        result = eliminateBD_resultASE(result);
        calculateXas(mXRay,result[0]);
        
    }
    public static void calculateXas(Double[][] matrix_XE,Double matrix_AS){
        System.out.print("\nX  | +a  | +s  |P(X,+a,+s)\n");
        System.out.printf("+x   | +a  | +s  |%f\n",matrix_XE[0][2]*matrix_AS);
        System.out.printf("-x   | +a  | +s  |%f\n",matrix_XE[1][2]*matrix_AS);
    }
    
    public static Double[] eliminateBD_resultASE(Double[] matrix_ASEBD){
        System.out.println("\nELIMINANDO B -> P(A,S,E,D)\n");
        Double[] matrix_ASE = new Double[2];
        
        matrix_ASE[0] = 0.;
        
        for(int i = 0; i < 4; i++){
            matrix_ASE[0] += matrix_ASEBD[i];
        }
        
        System.out.printf("%f\n",matrix_ASE[0]);
        
        return matrix_ASE;
    }
    
    public static Double[] calculateASEBD(Double[][] matrix_DEB,Double[] matrix_ASEB){
        System.out.println("\nCALCULANDO P(A,S,D,E,B)");
        Double[] matrix_ASEBD = new Double[4];
        
        for(int i = 0; i < 4; i++){
            if(i < 2){
                //System.out.printf("%f %f\n",matrix_DEB[i][3],matrix_ASEB[0]);
                matrix_ASEBD[i] = matrix_DEB[i][3] * matrix_ASEB[0];
            }else{
                //System.out.printf("%f %f\n",matrix_DEB[i][3],matrix_ASEB[1]);
                matrix_ASEBD[i] = matrix_DEB[i][3] * matrix_ASEB[1];
            }
            System.out.printf("%f\n",matrix_ASEBD[i]);
        }
        
        
        
        return matrix_ASEBD;
    }
    
    public static Double[] eliminateTL_resultASEB(Double[] matrix_ASETLB){
        System.out.println("\nELIMINANDO TL -> P(A,S,E,B)\n");
        Double[] matrix_ASEB = new Double[2];
    
        matrix_ASEB[0] = 0.;
        matrix_ASEB[1] = 0.;
        
        for(int i = 0; i < 8; i++){
            if(i % 2 == 0){
                matrix_ASEB[0] += matrix_ASETLB[i];
            }else{
                matrix_ASEB[1] += matrix_ASETLB[i];
            }
        }
        
        for(int i = 0; i < 2; i++){
            System.out.printf("%f\n",matrix_ASEB[i]);
        }
        
        return matrix_ASEB;
    }
    
    public static Double[] calculateASETLB(Double[][] matriz_ELT,Double[] matrix_ASTLB){
        System.out.println("\nCALCULANDO P(A,S,E,T,L,B)\n");
        Double[] matrix_ASETLB = new Double[8];
        
        for(int i = 0; i < 8; i++){
            matrix_ASETLB[i] = matriz_ELT[i/2][3] * matrix_ASTLB[i];
            System.out.printf("%f\n",matrix_ASETLB[i]);
        }
        
        return matrix_ASETLB;
    }
    
    public static Double[] calculateASTLB(Double[] matrix_TA,Double[] matrix_LBS){
        System.out.println("\nCALCULANDO P(A,S,T,L,B)\n");
        Double[] matrix_ASTLB = new Double[8];
        
        for(int i = 0; i < 8; i++){
            if(i < 4){
                matrix_ASTLB[i] = matrix_TA[0] * matrix_LBS[i];
            }else{
                matrix_ASTLB[i] = matrix_TA[1] * matrix_LBS[i-4];
            }
            System.out.printf("%f\n",matrix_ASTLB[i]);
        }
    
        return matrix_ASTLB;
    }
    
    public static Double[] evidenceS_resultLBs(Double[] matrix_LBS){
        Double[] matrix_LBs = new Double[4];
        
        for(int i = 0,j = 0; i < 8; i++){
            if(i % 2 == 0){
                matrix_LBs[j] = matrix_LBS[i];
                System.out.println(matrix_LBs[j]);
                j++;
            }        
        }   
        return matrix_LBs;
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
        
//        for(int i = 0; i < 8; i++){
//            for(int j = 0; j < 4; j++){
//                
//                System.out.print(matrix_LBS[i][j]+" ");
//            }
//            System.out.println("");
//        }
        
        return matrix_LBS;
    }
    
    public static Double[] calculateTA(Double[][] matrix_A, Double[][] matrix_T_A){
        System.out.println("\nCALCULANDO TA\n");
        Double[] matrix_TA = new Double[2];
        int index1,index2;
        
        for(int i = 0; i < 2; i++){
            if(i < 2){
                    index1 = 0;
            }else{
                index1 = 1;
            }
            matrix_TA[i] = matrix_A[index1][1] * matrix_T_A[i][2]; 
        }
               
        for(int i = 0; i < 2; i++){
            System.out.printf("%f \n",matrix_TA[i]);
        }
        
        return matrix_TA;
    }
}
