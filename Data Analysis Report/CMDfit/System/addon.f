        character*(*) function addon(ifname, ext)
                                                 
        character*(*) ifname, ext
                                 
                                 
       Integer Last, Dot, K, Len_1, Len_2, One
                                        
       Len_1 = Len( ifname )
       Len_2 = Len( Ext )
       Last = 0
       Dot = 0
c       Semi = 0
       Do K = 1, Len_1
         One = Ichar( ifname( K: K ) )
         If( One .eq. 47 ) Then
c                                   ! '/' End of Dir specification
           Dot = 0
         Else If( One .eq. 46 ) Then
c                                   ! '.' Extension if not Dir
           Dot = K
         Else If( One .gt. 32 ) Then
c                                   ! .not. ( ' ' End of file name )
           Last = K
         End If
       End Do
C
        addon=ifname
       if (Dot.eq.0 .and. ext(1:1).ne.' ' .and. ext(1:1).ne.char(0))
     &  addon = ifname( 1: Last ) // '.' // Ext
                                               
        end
