       blockdata blknam
                 
        character file*100, filex*3, lstfil*100
        logical inter
        integer lastf
                     
        common /arknmc/ file, filex, lstfil
        common /arknmn/ inter, lastf
                                                 
c       Initialise file name and extension.
        data file /'          '/
        data filex /'   '/
                          
c       Interactive mode is to be used by default.
        data inter /.true./
                           
c       The default file type to be written is ARK.
        data lastf /1/
                      
        end
