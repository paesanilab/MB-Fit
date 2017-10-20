"--Beginning of .vimrc file -- 

" This wraps lines at 80 characters 
 set textwidth=80

"set 20 lines to cursor
 set so=20

" creates a column showing 80 char length
" set colorcolumn=80

" allows you to see where you are in file
 set ruler

" " Sets tabs to be 2 characters instead of the default which is 8. 
 set tabstop=2
"
" " Number of spaces to use for each step of (auto)indent. 
 set shiftwidth=2
"
" " Tells vim to use spaces instead of tabs 
 set expandtab 
"
" " Tells vim to use c-style indenting 
 set cindent 
"
" " No mas press tab all day; the computers are gettings smarter
 set autoindent
"
"Line numbas
 set number
"

" Linebreak on 80 characters
 set lbr
 set tw=80
"
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
" => Colors and Fonts
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
" Enable syntax highlighting
syntax enable 

" highlights searched items
 set hlsearch

" makes searching better
 set incsearch

try
    colorscheme elflord
catch
endtry

"everything past column 80 is dark red
"let &colorcolumn=join(range(81,999),",")

"column at 80 
"let &colorcolumn="80,".join(range(120,999),",")

"highlight OverLength ctermbg=red ctermfg=white guibg=#59292
"match OverLength /\%81v.\+/

match Error /\%81v.\+/


"turn of error bells (sounds during mistakes)
set noerrorbells

" shows other files when using :e
set wildmenu
set wildmode=longest:full
set wildignore=*.o,*.bak,*.data,*.class

" Set 'comments' to format dashed lists in comments.
setlocal comments=sO:*\ -,mO:*\ \ ,exO:*/,s1:/*,mb:*,ex:*/,://
setlocal formatoptions+=ro



" Set extra options when running in GUI mode
"if has("gui_running")
"    set guioptions-=T
"    set guioptions-=e
"    set t_Co=256
"    set guitablabel=%M\ %t
"endif

" Set utf8 as standard encoding and en_US as the standard language
"set encoding=utf8

" Use Unix as the standard file type
"set ffs=unix,dos,mac
" -- End of .vimrc file --
