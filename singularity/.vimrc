inoremap jk <ESC>
syntax on
set number
set hlsearch
set ignorecase
set incsearch
set expandtab
set tabstop=4
set shiftwidth=4
set smarttab

augroup vimrc_todo
    au!
    au Syntax * syn match MyTodo /\v<(FIXME|NOTE[S]{0,1}|TODO|OPTIMIZE)/ containedin=ALL
augroup end
hi def link MyTodo Todo

if has("autocmd")
    augroup templates
        autocmd BufNewFile *.sh 0r ~/.vim/templates/skeleton.sh
        autocmd BufNewFile *.py 0r ~/.vim/templates/skeleton.py
"        autocmd BufNewFile *.R 0r ~/.vim/templates/skeleton.R
    augroup END
endif

map <F2> :retab <CR>

:nmap <c-s> :w<CR>
:imap <c-s> jk:w<CR>a
let @h = ':set ts=2 sts=2 noet | retab! | set ts=4 sts=4 et | retab!'
