#include<Rcpp.h>
using namespace Rcpp;
// Obtains the subset of values matching condition
// [[Rcpp::export]]


NumericVector STRPB_cpp(  NumericMatrix  X    ,
                              double        k    ,
                              int        BS   ,
                            double       rho
                          ) {
  Function  Uni("unique")                          ;
  Function  SAM("sample")                          ;
  Function  which("which")                         ;
  NumericMatrix           STR = (Uni(X))           ;   // unique values of X
  int                      n  =  X.nrow()          ;   // sample size
  int                      S  =  STR.nrow()        ;   // number of stratum
  int                    SIZE =  2*k*rho           ;   // size choose from sample form 1:BS
  NumericMatrix          IDSTR( n , S )            ;   // Create index for S stratum
  NumericVector           N_obs_S( S  )            ;   // number of sample in each stratum
  NumericVector           Treatment( n )           ;   // treatment
  NumericVector           b(BS)                    ;   // 1:BS
                          b =    seq_len(BS)-1     ;
                         
                         
  for( int i = 0 ; i < S ; i++ ){
   for( int j = 0 ; j < n ; j++  ){
     if(is_true(all(X(j,_)==STR(i,_) )) ){
       IDSTR ( j , i ) = 1                 ;
       } else{
       IDSTR ( j , i ) = 0                 ;
     }
    }
       N_obs_S[i] = sum( IDSTR(_,i) ) ;
 }
  
  NumericVector n_block =  ceiling( N_obs_S / BS )         ;
  
  for ( int i = 0 ; i < S ; i++ ){
    NumericVector    ID( N_obs_S[i] )                      ;
                     ID = which( IDSTR( _, i ) == 1 )      ;
                     ID = ID - 1                           ;
    for ( int j = 0 ; j < n_block[i] ; j++ ){
     NumericVector SEL =  SAM( b , Named("size") = SIZE )  ;
                   SEL =  SEL + j*BS                       ;
    for ( int t= 0 ; t< SIZE ; t++ )
      if( SEL[t] < ID.size() ){
         Treatment[ ID[SEL[t]]  ] = 1                      ;
        }
    }
  }

  return  Treatment ;
}

// [[Rcpp::export]]
NumericVector PS_cpp(  NumericMatrix  X ,
                       NumericVector  wt,
                       double        rho,
                       double         P    ) {
  int   n  = X.nrow();
  int   p  = X.ncol();
  NumericVector   Treatment = Rcpp::rbinom(  n , 1 ,  1/2 );
  NumericVector   Xcov(p)     ;
  NumericVector   Tmar(p)     ;
  double          delta       ;
  
  for( int i = 1 ; i < n ; i++ ){
    
                  Xcov  =  X.row(i)                        ;
    NumericMatrix  Xi   =  X(Rcpp::Range(seq(0,(i-1) )),_) ;
    NumericVector  Ti   =  Treatment[seq(0,(i-1))]         ;
    NumericVector  Dmar( p )                               ;
    
    for( int j = 0 ;j < p ; j++ ){
      NumericVector  mar_t     =  Ti[ Xi.column(j)==Xcov[j]] ;
                     Dmar[j]   =  sum(mar_t-rho)*wt[j]             ;
    }
                   delta =  sum(Dmar) ;

   if (delta>0)  {
     Treatment[i] = R::rbinom(1,1-P) ;
   } else if (delta<0){
     Treatment[i] = R::rbinom(1,P);
   } else{
     Treatment[i] = R::rbinom(1,1/2);
   }
  }
  return   Treatment ;
}


// [[Rcpp::export]]
NumericVector HH_cpp(  NumericMatrix  X ,
                       NumericVector  wt,
                       double        rho,
                       double         P     ) {
  int                     n = X.nrow();
  int                     p = X.ncol();
  NumericVector   Treatment = Rcpp::rbinom(  n , 1 ,  1/2 );
  NumericVector   Xcov(p)     ;
  NumericVector   Tmar(p)     ;
  double          Dstr        ;
  double          Dn          ;
  double          delta       ;
  
  for( int i = 1 ; i < n ; i++ ){
    
                   Xcov = X.row(i)                        ;
    NumericMatrix  Xi   = X(Rcpp::Range(seq(0,(i-1) )), _) ;
    NumericVector  Ti   = Treatment[seq(0,(i-1))]         ;
    NumericVector  Dmar( p )                              ;
    NumericVector  STR (i-1)                              ;
    
    for( int j = 0 ; j  < p    ; j++ ){
      NumericVector  mar_t     =  Ti[ Xi.column(j)==Xcov[j]] ;
                     Dmar[j]   =  sum( mar_t - rho )*wt[j+1]     ;
    }
    
    for( int t = 0 ; t < (i-1) ; t++ ){
      bool          De =  TRUE       ;
      for (int r = 0; r<p;r++){
        De = De & (Xi(t,r)==Xcov[r]);
      }
      
    if(De){
      STR[t]=t;
    }else{
      STR[t]=-1;
    }
     
    }
    
    NumericVector ST    = STR[STR>=0]                            ;
    NumericVector Tstr  = Ti[ST]                                 ;
                  Dstr  = sum(Tstr-rho) *wt[p+1]                 ;
                  Dn    = sum(Ti-rho)   *wt[0]                   ;
                  delta = Dn+ sum(Dmar)+Dstr                     ;
    
    if (delta>0)  {
      Treatment[i] = R::rbinom(1,1-P) ;
    } else if (delta<0){
      Treatment[i] = R::rbinom(1,P);
    } else{
      Treatment[i] = R::rbinom(1,1/2);
    }
  }
  return  Treatment ;
}





