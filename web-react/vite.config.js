import { defineConfig } from 'vite';
import react from '@vitejs/plugin-react';
export default defineConfig({
  plugins: [react()],
  base: './',
  build: {
    // The procedural 3D hero is isolated in a lazy chunk (~223 kB gzip).
    chunkSizeWarningLimit: 850,
  },
});
