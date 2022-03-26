       Subroutine RDAI_ERROR( Flag )
C
       Integer Flag
C
       If( Flag .eq. -1 ) Then
         Write( *, 101 )
101         Format( ' @ User request NOT to read a file' )
       Else If( Flag .eq. -2 ) Then
         Write( *, 102 )
102         Format( ' ERROR:: Input file not found in batch mode' )
       Else If( Flag .eq. -11 ) Then
         Write( *, 111 )
111         Format( ' ERROR:: Not an ARK file' )
       Else If( Flag .eq. -12 ) Then
         Write( *, 112 )
112         Format( ' ERROR:: Not a VAX (byte order) ARK file' )
       Else If( Flag .eq. -13 ) Then
         Write( *, 113 )
113         Format( ' ERROR:: Input FITS files is not one-dimensional' )
       Else If( Flag .eq. -14 ) Then
         Write( *, 114 )
114         Format
     &    (' ERROR:: Conflict in multiple definition of the data size')
       Else If( Flag .eq. -21 ) Then
         Write( *, 121 )
121         Format( ' ERROR:: Size of axis2 is not defined' )
       Else If( Flag .eq. -22 ) Then
         Write( *, 122 )
122         Format( ' ERROR:: Not a 1 dimensional file' )
       Else If( Flag .eq. -98 ) Then
         Write( *, 198 )
198         Format( ' ERROR:: This ASCII file too big' )
       Else If( Flag .eq. -99 ) Then
         Write( *, 199 )
199         Format( ' ERROR:: This ARK file is too big' )
       End If
C
       End
