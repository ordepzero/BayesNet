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
    }
    
    public static Double[] calculateASTLB(Double[] matriz_TA,Double[] matriz_LBS){
        System.out.println("\nCALCULANDO P(A,S,T,L,B)\n");
        Double[] matrix_ASTLB = new Double[8];
        
        for(int i = 0; i < 8; i++){
            if(i < 4){
                matrix_ASTLB[i] = matriz_TA[0] * matriz_LBS[i];
            }else{
                matrix_ASTLB[i] = matriz_TA[1] * matriz_LBS[i-4];
            }
            System.out.printf("%d %f\n",i,matrix_ASTLB[i]);
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
