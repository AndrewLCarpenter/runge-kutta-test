
  +  6   k820309    a          16.0        ���W                                                                                                           
       Newton.f90 NEWTON              NEWTON_ITERATION                      @                              
       WP                                                                                                       #         @                                                    	   #IPROB    #L    #EP    #DT    #NVECLEN    #TIME 	   #AI 
   #ICOUNT    #K                                                       
  @                                                    
  @                                                    
  @                                   
                
  @                                   
                
  @                                                    
  @                              	     
                
  @                              
     
                
D @                                                     D @                                           %         @   @                                                       #UVECITER    #K              
                                                    
              &                                                     
                                             (        `   @                                                              
    #BUILD_RNEWTON%UVEC    #EP    #DT    #TIME    #AI    #IPROB    #L    p          H r      7
S
r
   j            j                                      H r      7
S
r
   j            j                                                              @                                                
                &                                                     
  @                                   
                
                                      
                
  @                                   
                
  @                                   
                
  @                                                    
  @                                          #         @     @                                                #EP    #DT    #TIME    #AI    #IPROB    #L              
  @                                   
                
                                      
                
  @                                   
                
  @                                   
                
  @                                                    
  @                                          #         @     @                                                 #RNEWTON !          0  @                              !                   
               &                                           #         @     @                            "                 	   #IPROB #   #L $   #EP %   #DT &   #NVECLEN '   #TIME (   #AI )   #ICOUNT *   #K +             
  @                               #                     
  @                               $                     
  @                              %     
                
  @                              &     
                
                                  '                     
  @                              (     
                
  @                              )     
                
D                                 *                      
D                                 +            #         @     @                            ,                    #A -             
 @                              -                   
              &                   &                                           #         @     @                            .                    #MAT /   #NVECLEN 0   #RNEWTON 1             
D @                              /                   
               &                   &                                                     
  @                               0                     
D @                              1                   
               &                                           (        `   @                           2                                    
    #MAT 3     p        H r 4     7
S
l
O p          & &                 p          n                                          1  p         j            j                                  p          H r 4     7
S
l
O p          & &                 p          n                                          1  p         j            j                                    H r 4     7
S
l
O p          n                                      1  & &                 p          p         j            j                                      H r 4     7
S
l
O p          & &                 p          n                                          1  p         j            j                                    H r 4     7
S
l
O p          n                                      1  & &                 p          p         j            j                                                                        0  
                                3                   
              &                   &                                                         @                           4     SIZE               @                                SIZE    �         fn#fn    �   !   b   uapp(NEWTON    �   C   J  PRECISION_VARS "     p       WP+PRECISION_VARS !   �  �       NEWTON_ITERATION '   S  @   a   NEWTON_ITERATION%IPROB #   �  @   a   NEWTON_ITERATION%L $   �  @   a   NEWTON_ITERATION%EP $     @   a   NEWTON_ITERATION%DT )   S  @   a   NEWTON_ITERATION%NVECLEN &   �  @   a   NEWTON_ITERATION%TIME $   �  @   a   NEWTON_ITERATION%AI (     @   a   NEWTON_ITERATION%ICOUNT #   S  @   a   NEWTON_ITERATION%K    �  e       CHECK_EXIT $   �  �   a   CHECK_EXIT%UVECITER    �  @   a   CHECK_EXIT%K    �  �      BUILD_RNEWTON 5   d  �     BUILD_RNEWTON%UVEC+CONTROL_VARIABLES !   �  @   a   BUILD_RNEWTON%EP !   0  @   a   BUILD_RNEWTON%DT #   p  @   a   BUILD_RNEWTON%TIME !   �  @   a   BUILD_RNEWTON%AI $   �  @   a   BUILD_RNEWTON%IPROB     0	  @   a   BUILD_RNEWTON%L    p	  |       BUILD_JAC    �	  @   a   BUILD_JAC%EP    ,
  @   a   BUILD_JAC%DT    l
  @   a   BUILD_JAC%TIME    �
  @   a   BUILD_JAC%AI     �
  @   a   BUILD_JAC%IPROB    ,  @   a   BUILD_JAC%L    l  U       LU_SOLVER "   �  �   a   LU_SOLVER%RNEWTON !   M  �       NEWT_LINE_SEARCH '   �  @   a   NEWT_LINE_SEARCH%IPROB #   )  @   a   NEWT_LINE_SEARCH%L $   i  @   a   NEWT_LINE_SEARCH%EP $   �  @   a   NEWT_LINE_SEARCH%DT )   �  @   a   NEWT_LINE_SEARCH%NVECLEN &   )  @   a   NEWT_LINE_SEARCH%TIME $   i  @   a   NEWT_LINE_SEARCH%AI (   �  @   a   NEWT_LINE_SEARCH%ICOUNT #   �  @   a   NEWT_LINE_SEARCH%K    )  O       CONVERT_TO_CSR !   x  �   a   CONVERT_TO_CSR%A      k       QR_DECOMP    �  �   a   QR_DECOMP%MAT "   +  @   a   QR_DECOMP%NVECLEN "   k  �   a   QR_DECOMP%RNEWTON    �        MAT_INVERT      �   a   MAT_INVERT%MAT     �  =      MAT_INVERT%SIZE #   �  =      BUILD_RNEWTON%SIZE 