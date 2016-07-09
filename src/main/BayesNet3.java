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
        
        Double[] result = calculateLBS(mSmoke,mBronchitis,mCancer);
        result = normalize(result);
        Double[] result2 = calculateTA(mVisitAsia,mTuberculosis_VisitAsia);
        normalize(result2);
        result = calculateASTLB(result2,result);
        result = normalize(result);
        result = calculateSTLBE(mTBorCancer,result);
        result = eliminateTL_resultSEB(result);
        result = normalize(result);
        result = calculateSEBD(mDispnea_Bronchitis_TB_Cancer,result);
        result = eliminateEB_resultSD(result);
        result = normalize(result);
        calculateXas(mXRay,result);
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
    
    public static void calculateXas(Double[][] matrix_XE,Double[] matrix_AS){
        Double x1 = (matrix_XE[0][2]*matrix_AS[0])/(matrix_XE[0][2]*matrix_AS[0]+matrix_XE[1][2]*matrix_AS[1]);
        System.out.print("\nX  | +a  | +s  |P(X,+a,+s)\n");
        System.out.printf("+x   | +a  | +s  |%f\n",x1);
        System.out.printf("-x   | +a  | +s  |%f\n",1-x1);
    }
    
    public static Double[] eliminateEB_resultSD(Double[] matrix_SEBD){
        //System.out.println("\nELIMINANDO EB -> P(S,D)\n");
        Double[] matrix_SED = new Double[2];
        
        matrix_SED[0] = matrix_SEBD[0]+matrix_SEBD[2]+matrix_SEBD[4]+matrix_SEBD[6];
        matrix_SED[1] = matrix_SEBD[8]+matrix_SEBD[10]+matrix_SEBD[12]+matrix_SEBD[14];
        
//        for(int i = 0; i < 2;i++){
//            System.out.printf("%f\n",matrix_SED[i]);
//        }
        return matrix_SED;
    }
    
    public static Double[] calculateSEBD(Double[][] matrix_DEB,Double[] matrix_SEB){
        //System.out.println("\nCALCULANDO P(S,D,E,B)");
        Double[] matrix_SEBD = new Double[16];
        int i1 = 0;
        for(int i = 0; i < 16; i++){
            if(i == 0 || i == 8){
                i1 = 0;
            }else if(i == 1 || i == 9){
                i1 = 4;
            }else if(i == 2 || i == 10){
                i1 = 1;
            }else if(i == 3 || i == 11){
                i1 = 5;
            }else if(i == 4 || i == 12){
                i1 = 2;
            }else if(i == 5 || i == 13){
                i1 = 6;
            }else if(i == 6 || i == 14){
                i1 = 3;
            }else if(i == 7 || i == 15){
                i1 = 7;
            }
            matrix_SEBD[i] = matrix_DEB[i1][3] * matrix_SEB[i/2];
            //System.out.printf("%f\n",matrix_SEBD[i]);
        }
        
        return matrix_SEBD;
    }
    
    public static Double[] eliminateTL_resultSEB(Double[] matrix_STLBE){
        //System.out.println("\nELIMINANDO TL -> P(S,E,B)\n");
        Double[] matrix_SEB = new Double[8];
    
        matrix_SEB[0] = 0.;
        matrix_SEB[1] = 0.;
        matrix_SEB[2] = 0.;
        matrix_SEB[3] = 0.;
        matrix_SEB[4] = 0.;
        matrix_SEB[5] = 0.;
        matrix_SEB[6] = 0.;
        matrix_SEB[7] = 0.;
        
        for(int i = 0; i < 8; i++){
            if(i < 4){
                for(int j = i; j <= 12+i; j = j + 4){
                    matrix_SEB[i] = matrix_STLBE[j];
                }
            }else{
                for(int j = 12+i; j <= 24+i; j = j + 4){
                    matrix_SEB[i] = matrix_STLBE[j];
                }
            }
        }
        
//        for(int i = 0; i < 8; i++){
//            System.out.printf("%f\n",matrix_SEB[i]);
//        }
        
        return matrix_SEB;
    }
    
    public static Double[] calculateSTLBE(Double[][] matriz_ELT,Double[] matrix_ASTLB){
        //System.out.println("\nCALCULANDO P(S,E,T,L,B)\n");
        Double[] matrix_STLBE = new Double[32];
        int i,i1 = 0;
        
        for(int j = 0; j < 32; j++){
            i = j;
            if(j > 15){
                i = j - 16;
            }
            if(i == 0 || i == 2){
                i1 = 0;
            }else if(i == 1 || i == 3){
                i1 = 4;
            }else if(i == 4 || i == 6){
                i1 = 2;
            }else if(i == 5 || i == 7){
                i1 = 6;
            }else if(i == 8 || i == 10){
                i1 = 1;
            }else if(i == 9 || i == 11){
                i1 = 5;
            }else if(i == 12 || i == 14){
                i1 = 3;
            }else if(i == 13 || i == 15){
                i1 = 7;
            }
            if(j < 16){
                matrix_STLBE[j] = matrix_ASTLB[j] * matriz_ELT[i1][3];
                //System.out.printf("%d %f %f %f\n",i1,matrix_STLBE[j], matrix_ASTLB[j] , matriz_ELT[i1][3]);
            }else{
                matrix_STLBE[j] = matrix_ASTLB[j-16] * matriz_ELT[i1][3];
                //System.out.printf("%d %f %f %f\n",i1,matrix_STLBE[j],matrix_ASTLB[j-16], matriz_ELT[i1][3]);
            }
            //System.out.printf("%f\n",matrix_STLBE[j]);
        }
        
        return matrix_STLBE;
    }
    
    public static Double[] calculateASTLB(Double[] matrix_TA,Double[] matrix_LBS){
        //System.out.println("\nCALCULANDO P(A,S,T,L,B)\n");
        Double[] matrix_ASTLB = new Double[16];
        
        for(int i = 0; i < 16; i++){
            if(i < 8){
                matrix_ASTLB[i] = matrix_TA[0] * matrix_LBS[i];
            }else{
                matrix_ASTLB[i] = matrix_TA[1] * matrix_LBS[i-8];
            }
            //System.out.printf("%f\n",matrix_ASTLB[i]);
        }
    
        return matrix_ASTLB;
    }
    
    public static Double[] evidenceS_resultLBs(Double[] matrix_LBS){
        //System.out.println("\n\n");
        Double[] matrix_LBs = new Double[4];
        
        for(int i = 0,j = 0; i < 8; i++){
            if(i % 2 == 0){
                matrix_LBs[j] = matrix_LBS[i];
                //System.out.println(matrix_LBs[j]);
                j++;
            }        
        }   
        return matrix_LBs;
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
//            System.out.print(matrix_LBS[i]+" ");
//            System.out.println("");
//        }
        
        return matrix_LBS;
    }
    
    public static Double[] calculateTA(Double[][] matrix_A, Double[][] matrix_T_A){
        //System.out.println("\nCALCULANDO TA\n");
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
               
//        for(int i = 0; i < 2; i++){
//            System.out.printf("%f \n",matrix_TA[i]);
//        }
        
        return matrix_TA;
    }
}
