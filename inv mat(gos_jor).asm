; Inversion Of Matrix nxn 
; by ShahroOz

; This code uses Gauss Jordan algorithm
; for calculation

; HELP:
;    place your matrix in mat and the size
;    of that in N. just run it. 


data segment
    temp1 dw ?
    temp2 dw ?
    temp3 dw ?
    l dw ?
    numstr db '$$$$$'
    manfi dw -1
    N dw 3
    ;mat dw 1,2,0,0,-3,4,0,0
    ;mat dw 4,1,2,12
    mat dw -5,-1,8,1,1,7,-1,10,-1
    new_mat dw 100 dup(?)
data ends

stack segment
    mystack db 100 dup(?)
stack ends

myprog segment
    main proc far
        assume cs: myprog, ds: data, ss: stack
        mov ax, data
        mov ds, ax
        lea sp, mystack+100
        
        lea si, mat
        mov bp, n
        lea di, new_mat
        call make
        lea si, new_mat        
        call func1; setting adjacent matrix
        call func2; partial pivoting
        call gos_jor; do the magic... :))
        
        
        call print
        
        
        hlt
        
    main endp
    
    make proc near
        
        push bp
        push di
        push si
        push ax
        push bx
        push cx
        
        xor ax, ax;as i
  loopm1:cmp ax, N
        jge outm1
        xor bx,bx
  loopm2: cmp bx, N
        jge outm2
        mov cx, [si]
        mov [di], cx
        inc bx
        add si, 2
        add di, 2
        jmp loopm2
        
        
   outm2:
        xor bx, bx
  loopm3: cmp bx, N
        jge outm3
        mov cx, 0
        mov [di], cx
        inc bx
        ;add si, 2
        add di, 2
        jmp loopm3
        
        
   outm3: inc ax
        jmp loopm1
        
   outm1:
        
        pop cx
        pop bx
        pop ax
        pop si
        pop di
        pop bp
        
        ret
        
    make endp
    
    print proc near
        
        push bp
        push ax
        
        mov ax, 2
        mul bp
        mov bp, ax
        
        xor ax, ax
  loopp1: cmp ax, n
        jge outp1
        mov bx, n
  loopp2: cmp bx, bp
        jge outp2
        push bx
        
        push ax
        push dx
        push si
        mov ax, ax
        mov dx, bx
        mov si, si
        call aij
        mov ax, [si]
        mov l, ax
        pop si
        pop dx
        pop ax
        
        
        push ax
        push dx
        push si
        cmp l, 0
        jge ok1
        mov ah, 2
        mov dx, '-'
        int 21h
        push ax
        push dx
        mov ax, l
        xor dx,dx
        imul manfi
        mov l, ax
        pop dx
        pop ax
   ok1: push cx
        push si
        mov ax,l
        lea si, numstr
        call s2n
        pop si
        pop cx
        
        pop si
        pop dx
        pop ax
        
        push ax
        push dx
        mov ah, 2
        mov dx, '/'
        int 21h
        pop dx
        pop ax
        
        
        push ax
        push dx
        push si
        mov ax, ax
        mov dx, ax
        mov si, si
        call aij
        mov ax, [si]
        mov l, ax
        pop si
        pop dx
        pop ax
        
        push ax
        push dx
        push si
        cmp l, 0
        jge ok2
        mov ah, 2
        mov dx, '-'
        int 21h
        push ax
        push dx
        mov ax, l
        xor dx,dx
        imul manfi
        mov l, ax
        pop dx
        pop ax
   ok2: push cx
        push si
        mov ax,l
        lea si, numstr
        call s2n
        pop si
        pop cx
        pop si
        pop dx
        pop ax
         
        push ax
        push dx
        mov ah, 2
        mov dx, ' '
        int 21h
        pop dx
        pop ax
        
        pop bx
        
        inc bx
        jmp loopp2
        
        
 outp2: push ax
        push dx
        mov dx,13
        mov ah,2
        int 21h  
        mov dx,10
        mov ah,2
        int 21h
        int 21h
        pop dx
        pop ax
        
        inc ax
        jmp loopp1
    
 outp1: pop ax
        pop bp
        ret 
        
    print endp
    
    
    func2 proc near
        
        push bp
        push di
        push ax
        push dx
        push si
        
        dec bp
        mov ax, bp
        inc bp
        
        push ax
        mov ax, 2
        mul bp
        mov bp, ax
        pop ax
        
   loop14:
        cmp ax, 0
        jle out14
        
        push si
        push ax
        push dx
        mov ax, ax
        dec ax
        mov dx, 0
        call aij
        mov di, word ptr [si]
        mov temp1, di
        pop dx
        pop ax
        pop si 
        
        push si
        push ax
        push dx
        mov ax, ax
        mov dx, 0
        call aij
        mov di, word ptr [si]
        mov temp2, di
        pop dx
        pop ax
        pop si
        
        mov di, temp1
        cmp di, temp2
        jge out15 
        xor dx, dx 
   loop16: cmp dx, bp
        jge out16
        push si
        push ax
        push dx
        mov ax, ax
        mov dx, dx
        call aij
        mov di, word ptr [si]
        mov temp1, di
        pop dx
        pop ax
        pop si
        
        push si
        push ax
        push dx
        mov ax, ax
        dec ax
        mov dx, dx
        call aij
        mov di, word ptr [si]
        mov temp2, di
        pop dx
        pop ax
        pop si
               
        push si
        push ax
        push dx
        mov ax, ax
        mov dx, dx
        call aij         
        mov di, temp2 
        mov word ptr [si], di
        pop dx
        pop ax
        pop si  
        
        push si
        push ax
        push dx
        mov ax, ax
        dec ax
        mov dx, dx
        call aij         
        mov di, temp1 
        mov word ptr [si], di
        pop dx
        pop ax
        pop si      
        
        inc dx
        jmp loop16
        
        
    out16:
        
    out15: dec ax
        jmp loop14   
        
    out14:
     
        pop si
        pop dx
        pop ax
        pop di
        pop bp
        
        ret
        
    func2 endp
    
    func1 proc near;si as mat, bp as n
        
        push ax
        push dx
        push di
        push si
        push bp
        
        mov ax, 2
        xor dx, dx
        mul bp
        mov bp, ax;bp=2*n
        
        mov ax, 0
  loop21: cmp ax, n
        jge out21
        mov dx, 0
  loop22: cmp dx, bp
        jge out22
        mov di, ax
        add di, n
        cmp di, dx
        jne out23
        
        push si
        push ax
        push dx
        mov ax, ax
        mov dx, dx
        call aij
        mov word ptr [si], 1
        pop dx
        pop ax
        pop si
        
   out23: inc dx
        jmp loop22    
        
   out22: inc ax
        jmp loop21    
        
   out21: pop bp
        pop si     
        pop di
        pop dx
        pop ax
        ret         
          
    func1 endp
    
    aij proc near; si as arr, ax as i, dx as j
              
        add si, dx
        add si, dx
        xor dx, dx
        mul N
        add si, ax
        add si, ax
        add si, ax
        add si, ax
        
        ret
        
    aij endp
        
    gos_jor proc near; si as arr, bp as var
        
        push dx
        push bx
        push di
        push si
        push bp
        
        xor dx, dx
        mov ax, 2
        mul bp
        mov bp, ax;bp=2*n
        
        
        xor dx, dx
 loop1: cmp dx, N
        jge out1
        xor bx, bx
 loop2: cmp bx, bp
        jge out2
        
        push ax
        push dx
        push si
        mov ax, bx
        mov dx, dx
        mov si, si
        call aij
        mov ax, [si]
        mov l, ax
        pop si
        pop dx
        pop ax
        
        xor di, di
 loop3: cmp di, bp
        jge out3
        cmp bx, dx
        je out4
        push ax
        push dx
        push si
        mov ax, dx
        mov dx, dx
        mov si, si
        call aij
        mov ax, [si]
        mov temp1, ax
        pop si
        pop dx
        pop ax
        
        push ax
        push dx
        push si
        mov ax, bx
        mov dx, di
        mov si, si
        call aij
        mov ax, [si]
        pop si
        xor dx, dx
        imul temp1
        mov temp2, ax
        pop dx
        pop ax
        
        push ax
        push dx
        push si
        mov ax, dx
        mov dx, di
        mov si, si
        call aij
        mov ax, [si]
        pop si
        xor dx, dx
        imul l
        mov temp3, ax
        pop dx
        pop ax
        
        push ax
        mov ax, temp2
        sub ax, temp3
        mov temp3, ax
        push ax
        push dx
        push si
        mov ax, bx
        mov dx, di
        mov si, si
        call aij
        mov ax, temp3
        mov [si], ax
        pop si
        pop dx
        pop ax
        pop ax
        
        
  out4: inc di
        jmp loop3     
        
  out3: inc bx
        jmp loop2     
        
  out2: inc dx
        jmp loop1
        
  out1: pop bp
        pop si
        pop di
        pop bx
        pop dx
        
        ret    
    
    gos_jor endp
    
    s2n proc near
        call dollars 
        mov bx, 10  
        mov cx, 0   
loops11:       
        mov dx, 0   
        div bx      
        push dx      
        inc cx      
        cmp ax, 0   
        jne loops11   
loops21:  
        pop dx        
        add dl, 48  
        mov [si], dl
        inc si
        loop loops21
  
        mov ah, 9
        mov dx, offset numstr
        int 21h   

        ret
    s2n endp    
    
dollars proc near                 
        mov cx, 5
        mov di, offset numstr
loopd:  mov bl, '$'
        mov [ di ], bl
        inc di
        loop loopd

  ret
dollars endp
    
    
myprog ends

    end main
    