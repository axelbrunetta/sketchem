GRADIENT_BUTTON_CSS = """
    
                button:not([disabled]) {
    
                    background: linear-gradient(90deg, #0066cc, #4da6ff, #0066cc) !important;
                    background-size: 200% 100% !important;
                    color: white !important;
                    padding: 15px 32px !important;
                    font-size: 16px !important;
                    font-weight: bold !important;
                    border-radius: 12px !important;
                    border: none !important;

                    
                    transition: background-position 0.5s ease,
                                transform 0.3s ease, 
                                box-shadow 0.3s ease !important;
                }


                button:is(:hover, :focus-visible):not([disabled]) {
                    
                    background-position: 100% 0 !important;

                    
                    
                    transform: scale(1.05) !important; 
                    
                    box-shadow: 0 0 20px 5px rgba(77, 166, 255, 0.6) !important;
                

                }

                button[disabled] {
                
                    background-color: #cccccc !important;
                    color: #666666 !important;
                    cursor: not-allowed !important;
                }
                
                """