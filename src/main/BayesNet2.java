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
public class BayesNet2 {
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
        
        Double[][] result = calculateTA(mVisitAsia,mTuberculosis_VisitAsia);
        
        System.out.println("T  |  P(T)");
        System.out.println("+t |  "+(result[0][2]+result[2][2]));
        System.out.println("-t |  "+(result[1][2]+result[3][2]));
    }
    
    
    public static Double[][] calculateTA(Double[][] matrix_A, Double[][] matrix_T_A){
        //System.out.println("\nCALCULANDO TA\n");
        Double[][] matrix_TA = new Double[4][3];
        int index1,index2;
        
        for(int i = 0; i < 4; i++){
            if(i < 2){
                matrix_TA[i][0] = 1.;
            }else{
                matrix_TA[i][0] = 0.;
            }
        }
        
        for(int i = 0; i < 4; i++){
            if(i % 2 == 1){
                matrix_TA[i][1] = 0.;
            }else{
                matrix_TA[i][1] = 1.;
            }
        }
        
        for(int i = 0; i < 4; i++){
            if(i < 2){
                    index1 = 0;
            }else{
                index1 = 1;
            }
            matrix_TA[i][2] = matrix_A[index1][1] * matrix_T_A[i][2]; 
        }
               
//        for(int i = 0; i < 4; i++){
//            for(int j = 0; j < 3; j++){
//                
//                System.out.print(matrix_TA[i][j]+" ");
//            }
//            System.out.println("");
//        }
        
        return matrix_TA;
    }
}
