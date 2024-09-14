from flask import Flask, send_from_directory

from chemquest_website import app

@app.route('/get-audio/<filename>')
def get_audio(filename):
    audio_directory = 'static/audio'
    return send_from_directory(audio_directory, filename)
