/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package main;

import entities.Table;

/**
 *
 * @author PeDeNRiQue
 */
public class BayesNet {
    
    static String[] hVistitAsia = {"A","P(A)"};
    static Double[][] mVisitAsia = {{1., 0.01},
                                    {0., 0.99}};
    
    static String[] hTuberculosis_VistiAsia = {"A","T","P(T|A)"};
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
        // TODO code application logic here
        
        
        Table nVisitAsia = new Table(hVistitAsia,mVisitAsia);
        Table nTuberculosis_VisitAsia = new Table(hTuberculosis_VistiAsia,mTuberculosis_VisitAsia);
     
        joinTable(nVisitAsia,nTuberculosis_VisitAsia);
    }
    
    public static Double[][] joinTable(Table node1, Table node2){
        Double[][] resultTable;
        Integer[][] indexes = {{-1,-1,-1,-1,-1,-1},{-1,-1,-1,-1,-1}};
        String[] headers = new String[5];
        
        Double[][] m1 = node1.getMatrix();
        Double[][] m2 = node2.getMatrix();
        
        String[] h1 = node1.getHeaders();
        String[] h2 = node2.getHeaders();
        int index;
        
        for(int i = 0; i < node1.size(); i++){
            index = node2.hasHeader(h1[i]);
            if(index != -1){
                indexes[0][i] = i;
                indexes[1][i] = index;
            }
        }
        
        for(int i = 0; i < 5; i++){
            System.out.println(indexes[0][i]+" "+indexes[1][i]+"-"+headers[i]+"-");
        }
        
        int hasHeader;
        int headersSize = 0;
        for(int i = 0; i < h1.length - 1; i++){
            headers = hasHeader(headers,h1[i]);
        }
        
        for(int i = 0; i < h2.length - 1; i++){
            headers = hasHeader(headers,h2[i]);
        }
        
        for(int i = 0; i < headers.length; i++){
            System.out.println(headers[i]);
            if(headers[i] == null){
                headersSize = i;
            }
        }
        
        //CALCULA O TAMANHO DA TABELA RESULTANTE BASEADO NO NÚMERO DE COLUNAS (HEADERSSIZE)
        //O "+1" É PRA ADICIONAR A COLUNA DA PROBABILIDADE
        //E DE LINHAS 2^n
        resultTable = new Double[(int) Math.pow(2,headersSize)][headersSize+1];
        
        return resultTable;
    }
    
    public static String[] hasHeader(String[] headers,String header){
        for(int i = 0; i < headers.length; i++){
            if(headers[i] == null){
                headers[i] = header;
                return headers;
            }else if(headers[i].equals(header)){
                return headers;
            }
        }
        
        return headers;
    }
    

}
