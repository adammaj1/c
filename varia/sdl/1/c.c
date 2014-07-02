/*
https://en.wikibooks.org/wiki/X_Window_Programming/SDL

*/

// Headers
#include "<SDL/SDL.h>"
 
// Main function
int main( int argc, char* argv[] )
{
    // Initialize SDL
    if( SDL_Init( SDL_INIT_EVERYTHING ) == -1 )
        return 1;
 
    // Delay 2 seconds
    SDL_Delay( 2000 );
 
    // Quit SDL
    SDL_Quit();
 
    // Return
    return 0;
}
